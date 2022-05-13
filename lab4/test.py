
a = (5,5)
b = (0,0)
c = (3,2)

def cross():
    tmp = (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    if tmp==0:
        print('on line')
    elif tmp>0:
        print('left')
    else:
        print('right')

cross()