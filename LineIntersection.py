class Point():
    def __init__(self,x,y):
        self.x=x
        self.y=y

def OnSegment(p,q,r):
    if q.x<=max(p.x,r.x) and q.x>=min(p.x,r.x) and q.y<=max(p.y,r.y) and q.y>=min(p.y,r.y) :
        return True
    return False

def Orientation(p,q,r):
    val=(q.y-p.y)*(r.x-q.x)-(q.x-p.x)*(r.y-q.y)
    if val==0:
        return 0
    elif val>0:
        return 1 
    else:
        return 2
def IfLineIntersect(edge1,edge2):

    p1=Point(edge1[0],edge1[1])
    q1=Point(edge1[2],edge1[3])
    p2=Point(edge2[0],edge2[1])
    q2=Point(edge2[2],edge2[3])


    o1=Orientation(p1,q1,p2)
    o2=Orientation(p1,q1,q2)
    o3=Orientation(p2,q2,p1)
    o4=Orientation(p2,q2,q1)

    if (o1!=o2) and (o3!=o4):
        return True

    if (o1==0) and OnSegment(p1,p2,q1):
        return True
    if (o2==0) and OnSegment(p1,q2,q1):
        return True
    if (o3==0) and OnSegment(p2,p1,q2):
        return True
    if (o4==0) and OnSegment(p2,q1,q2):
        return True
    return False

# print IfLineIntersect((0,0,0,10),(0,0,0,-3))
     