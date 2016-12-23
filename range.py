y0 = 0.534327011356
d0 = 0.0382009329056
y1 = 0.553340049549
d1 = 0.0102917149077
y2 = 0.545450295339
d2 = 0.0103656866711
y3 = 0.596134682227
d3 = 0.0988576842745
y4 = 0.523702787918
d4 = 0.02957239783
y5 = 0.536227894832
d5 = 0.0349878620091

s0 = 0.00796544514125

ys = [y0,y1,y2,y3,y4,y5]
ds = [d0,d1,d2,d3,d4,d5]

maxes = []
mins  = []

for i in range(len(ys)):
    nmx = ys[i] + ds[i]
    nmn = ys[i] - ds[i]
    maxes.append(nmx)
    mins.append(nmn)

fmx = max(maxes)
fmn = min(mins)

print(maxes)
print(mins)

print('The final interval is: ',(fmx,fmn))
print('---------------------------------')
print((y0+d0+s0),(y0-d0-s0))
print((y0+d0+(2.*s0)),(y0-d0-(2.*s0)))
print((y0+d0+(3.*s0)),(y0-d0-(3.*s0)))
print('---------------------------------')
