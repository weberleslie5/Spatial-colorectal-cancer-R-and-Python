from PIL import Image, ImageChops
def trim(im):
    bg = Image.new("RGB", im.size, "white")
    return im.crop(ImageChops.difference(im.convert("RGB"), bg).getbbox())
base = "/home/weberl/projects/crc_visiumhd/figures/"
imgs = [trim(Image.open(base + f"rvspy_{s}.png")) for s in ["P1CRC", "P2CRC", "P5CRC"]]
w = min(i.width for i in imgs)
imgs = [i.resize((w, int(i.height * w / i.width))) for i in imgs]
pad = 40
out = Image.new("RGB", (w + pad * 2, sum(i.height for i in imgs) + pad * 4), "white")
y = pad
for i in imgs:
    out.paste(i, (pad, y)); y += i.height + pad
out.save(base + "fig_rvspy.png")
print("stitched")
