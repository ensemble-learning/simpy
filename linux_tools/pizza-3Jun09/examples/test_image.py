# simple test of image tool
# requires files/bucky*.png

i = image("files/bucky*.png")
i.convert("files/bucky*.png","tmp*.gif")
i.montage("","files/bucky*.png","tmp*.gif","tmpnew*.png")
i.view("")

print "all done ... type CTRL-D to exit Pizza.py"
