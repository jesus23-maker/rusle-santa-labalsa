"""Decode the exported RUSLE class and hotspot layers into arrays (results/grid.npz)."""
import numpy as np
from PIL import Image
cls_img = np.array(Image.open('data/layers/erosion_RUSLE.png').convert('RGBA'))
hot = np.array(Image.open('data/layers/hotspots_P90.png').convert('RGBA'))
COLS = np.array([[0x8c,0x6d,0x31],[0xb1,0x5c,0xe6],[0xf4,0xa5,0x82],[0x4d,0xaf,0x4a],[0xff,0x9f,0x1c],[0xd7,0x19,0x1c]])
valid = cls_img[..., 3] > 0
d = ((cls_img[..., None, :3].astype(int) - COLS[None, None]) ** 2).sum(-1)
cls = d.argmin(-1)
np.savez('results/grid.npz', cls=cls, valid=valid, hot=(hot[..., 3] > 0) & valid)
print('valid pixels', valid.sum(), '| basin hotspot fraction', round(((hot[..., 3] > 0) & valid).sum() / valid.sum(), 4))
