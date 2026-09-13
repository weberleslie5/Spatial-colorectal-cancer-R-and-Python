from PIL import Image, ImageChops
def trim(im):                                               # crop white borders
    bg = Image.new("RGB", im.size, "white")
    box = ImageChops.difference(im.convert("RGB"), bg).getbbox()
    return im.crop(box)
base = "/home/weberl/projects/crc_visiumhd/figures/"
names = ["joint_panel_P1CRC.png", "joint_panel_P2CRC.png", "joint_panel_P5CRC.png", "joint_legend.png"]
imgs = [trim(Image.open(base + n)) for n in names]
h = min(i.height for i in imgs[:3])                         # common height from the panels
imgs = [i.resize((int(i.width * h / i.height), h)) for i in imgs]
pad = 40
out = Image.new("RGB", (sum(i.width for i in imgs) + pad * 5, h + pad * 2), "white")
x = pad
for i in imgs:
    out.paste(i, (x, pad)); x += i.width + pad
out.save(base + "fig3a_joint.png")
print("stitched")
