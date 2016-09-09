import pygal
import hashlib
import random


def removeTop(mA):
    if mA.startswith("<?xml"):
        mA = mA[mA.index('?>') + 2:]
    return mA

bar = pygal.Line()
bar.add('A', xrange(10))

sq = pygal.Line()
sq.add('B', [i**2 for i in xrange(10)])

sq.render_to_file('re.svg')

f = "res.svg"
r = ""
u = ""
i = [[bar, sq]]
xD = 100 / len(i)
yD = 100 / len(i[0])
for x, vx in enumerate(i):
    for y, vy in enumerate(vx):
        d = hashlib.sha256(str(random.getrandbits(10))).hexdigest()
        r += '<svg id="{0}">{1}</svg>'.format(d, removeTop(vy.render()))
        u += '<use xlink="#{0}" width="{1}%" height="{2}%" x="{3}%" y="{4}%"/>'.format(d,xD,yD,x*xD,y*yD)
with open(f,'w') as fd:
    fd.write('<?xml version=\'1.0\' encoding=\'utf-8\'?>\n<svg xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg"><defs>{0}</defs>{1}</svg>'.format(r,u))
