import ee, math, json, time
ee.Initialize(project='rusle-santa-unalm')
OUT = {}
merit = ee.Image('MERIT/Hydro/v1_0_1')
def basin_from(lon, lat):
    pt = ee.Geometry.Point([lon, lat])
    s = merit.select('upa').addBands(ee.Image.pixelLonLat()).reduceRegion(ee.Reducer.max(3), pt.buffer(3000), 90).getInfo()
    sp = ee.Geometry.Point([s['max1'], s['max2']])
    hb = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_12')
    reg = hb.filterBounds(ee.Geometry.Rectangle([-78.6, -10.4, -76.9, -8.2]))
    pairs = reg.reduceColumns(ee.Reducer.toList(2), ['HYBAS_ID', 'NEXT_DOWN']).get('list').getInfo()
    ch = {}
    for h, n in pairs: ch.setdefault(int(n), []).append(int(h))
    st = [int(hb.filterBounds(sp).first().get('HYBAS_ID').getInfo())]; ids = []
    while st:
        h = st.pop(); ids.append(h); st.extend(ch.get(h, []))
    b = reg.filter(ee.Filter.inList('HYBAS_ID', ids)).union(30).geometry()
    return b, s, len(ids)
basin, s0, n0 = basin_from(-77.83, -8.86)
OUT['paper_outlet'] = dict(snap=[s0['max1'], s0['max2']], upa=s0['max'], n=n0, area=basin.area(30).divide(1e6).getInfo())
# Egenor (1999) coordinates of La Balsa: 77 49 47 W, 8 52 27 S
lb = ee.Geometry.Point([-(77 + 49/60 + 47/3600), -(8 + 52/60 + 27/3600)])
OUT['upa_at_egenor_point_max500m'] = merit.select('upa').reduceRegion(ee.Reducer.max(), lb.buffer(500), 90).getInfo()
dem = ee.Image('NASA/NASADEM_HGT/001').select('elevation')
OUT['dem_at_egenor_point'] = dem.reduceRegion(ee.Reducer.mean(), lb.buffer(100), 30).getInfo()
CRS, SCALE = 'EPSG:32718', 30
chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD').select('precipitation')
def R4_months(d0, d1, nyears):
    def f(m):
        m = ee.Number(m)
        return chirps.filterDate(d0, d1).filter(ee.Filter.calendarRange(m, m, 'month')).sum().divide(nyears).rename('P')
    Pm = ee.ImageCollection(ee.List.sequence(1, 12).map(f))
    Pa = Pm.sum()
    mfi = Pm.map(lambda im: im.pow(2)).sum().divide(Pa)
    return mfi.pow(0.552).multiply(161).rename('R4').resample('bilinear'), Pa
def sg(v): return ee.Image(f'projects/soilgrids-isric/{v}_mean').select(f'{v}_0-5cm_mean')
SAN = sg('sand').divide(10); SIL = sg('silt').divide(10); CLA = sg('clay').divide(10); OC = sg('soc').divide(100)
SN1 = SAN.divide(100).multiply(-1).add(1)
w1 = SAN.multiply(-0.256).multiply(SIL.divide(100).multiply(-1).add(1)).exp().multiply(0.3).add(0.2)
w2 = SIL.divide(CLA.add(SIL)).pow(0.3)
w3 = OC.multiply(0.25).divide(OC.add(OC.multiply(-2.95).add(3.72).exp())).multiply(-1).add(1)
w4 = SN1.multiply(0.7).divide(SN1.add(SN1.multiply(22.9).add(-5.51).exp())).multiply(-1).add(1)
K = w1.multiply(w2).multiply(w3).multiply(w4).multiply(0.1317)
th = ee.Terrain.slope(dem).multiply(math.pi / 180); sinth = th.sin(); s_pct = th.tan().multiply(100)
beta = sinth.divide(0.0896).divide(sinth.pow(0.8).multiply(3).add(0.56)); m_r = beta.divide(beta.add(1))
LS = ee.Image(30 / 22.13).pow(m_r).multiply(sinth.multiply(10.8).add(0.03).where(s_pct.gte(9), sinth.multiply(16.8).subtract(0.5)))
def s2p(im):
    scl = im.select('SCL'); ok = scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    return ee.Image(im.select(['B3', 'B4', 'B8', 'B11']).divide(10000).updateMask(ok).copyProperties(im, ['system:time_start']))
s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(basin).filterDate('2023-01-01', '2024-12-31').filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40)).map(s2p)
s2a = s2.median(); s2d = s2.filter(ee.Filter.calendarRange(6, 8, 'month')).median()
NDVI = s2a.normalizedDifference(['B8', 'B4']); NDSI = s2d.normalizedDifference(['B3', 'B11'])
C = NDVI.expression('exp(-2*n/(1-n))', {'n': NDVI.max(0).min(0.99)}).min(1)
wc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map')
glac = NDSI.gte(0.4).And(s2d.select('B4').gt(0.2)); water = wc.eq(80).Or(NDVI.lt(0).And(NDSI.lt(0.4)))
mask = glac.Or(water).Or(wc.eq(70)).Not()
KLSC = K.multiply(LS).multiply(C).updateMask(mask).clip(basin)
def meanA(R):
    return R.multiply(KLSC).reduceRegion(ee.Reducer.mean(), basin, SCALE, crs=CRS, maxPixels=1e13, tileScale=8).getInfo()
def meanP(Pa):
    return Pa.reduceRegion(ee.Reducer.mean(), basin, 5000).getInfo()
Rc, Pc = R4_months('1995-01-01', '2025-01-01', 30)
OUT['clim'] = dict(A=meanA(Rc), P=meanP(Pc))
print('clim', OUT['clim'], flush=True)
for start_m, tag in [(9, 'SepAug'), (10, 'OctSep')]:
    for y in range(2002, 2009):
        d0 = f'{y}-{start_m:02d}-01'; d1 = f'{y+1}-{start_m:02d}-01'
        R, Pa = R4_months(d0, d1, 1)
        OUT[f'{tag}_{y}-{y+1}'] = dict(A=meanA(R), P=meanP(Pa))
        print(tag, y, OUT[f'{tag}_{y}-{y+1}'], flush=True)
print('RESULT_JSON ' + json.dumps(OUT))
