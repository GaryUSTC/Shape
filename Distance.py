import math
import numpy as np
# import skimage
# from skimage import draw
import DPP 
import LineIntersection
import copy 
Width=600
Height=400
import time
from sympy import *
# import sympy


def GridentDescent(center_vectors,polygons,learning_rate,value,Groups):

  
    t=0.5
    flag=False

    current_vectors=copy.copy(center_vectors)
    
    
    for i in range(0,len(center_vectors)):
        
        c1=copy.copy(center_vectors)
        c1[i]+=t

        '''caculate the derivative, then use gradient descent]
            here we use approximate value instead of the accurate derivative'''
        det=(ObjectiveFunction(tuple(c1),polygons,Groups)-ObjectiveFunction(tuple(center_vectors),polygons,Groups))/t 


        center_vectors[i]+=(learning_rate[i]*det)

        # if det>0:
        #     # learning_rate[i]*=1.5
        #     center_vectors[i]+=(learning_rate[i]*det)                
        # else:
        #     # learning_rate[i]*=0.8
        #     center_vectors[i]+=(learning_rate[i]*det) 


        
        '''if polygon is going to overlap with others, reduce the learning_rate, i.e, iteration pace'''
        for j in range(0,len(polygons)):
            if i/2!=j and Groups[i/2]==Groups[j]:
                if isPolygonsOverlap1(center_vectors[i/2*2],center_vectors[i/2*2+1],center_vectors[j*2],center_vectors[j*2+1],polygons[i/2],polygons[j]):
                    center_vectors[i]=current_vectors[i]+(learning_rate[i]/10*det)    
                    break

        ''' if polygon is out of bound, recall this iteration'''
        if isPolygonOutOfBound(center_vectors[i/2*2],center_vectors[i/2*2+1],polygons[i/2]):
            center_vectors[i]=current_vectors[i]

    return ObjectiveFunction(tuple(center_vectors),polygons,Groups)

def isPolygonsOverlap1(x1,y1,x2,y2,p1,p2):
    '''to check whether two polygons (p1,p2) with translations(x1,y1,x2,y2) overlap'''
    temp1=[]
    for i in p1:
        if i%2==0:
            temp1.append(i+x1)
        else:
            temp1.append(i+y1)
    temp2=[]
    for i in p2:
        if i%2==0:
            temp2.append(i+x2)
        else:
            temp2.append(i+y2)
    return isPolygonsOverlap(temp1,temp2)



def isPolygonOutOfBound(x,y,polygon):
    ''' to check whether a polygon with translation(x,y) is out of bound'''
    for i in range(0,len(polygon)):
        if i%2==0:
            a=polygon[i]+x
            if a<0 or a>Width:
                return True
        else:
            a=polygon[i]+y
            if a<0 or a>Height:
                return True
    return False


def ObjectiveFunction(center_vectors,polygons,Groups):
    ''' this function is in order to make translate'''
    new_polygons=[]
    for i in range(0,len(center_vectors)/2):
        temp=[]
        for j in range(0,len(polygons[i])):
            if(j%2==0):
                temp.append(polygons[i][j]+center_vectors[2*i])
            else:
                temp.append(polygons[i][j]+center_vectors[2*i+1])
        new_polygons.append(tuple(temp))
    
    # return SumOfDistances(new_polygons)+0.5*MinOfDistances(new_polygons)

   
    result = test(new_polygons,Groups)
    return result



def test(polygons,Groups):
    '''the true objective function
    for each polygon, we caculate the shortest distance of different group subtract the longest distance of the same group'''
    sa=0.0
    di=0.0
    for i in range(0,len(polygons)):
        dis_same=[]
        dis_diff=[]
        temp=[]
        for j in range(0,len(polygons)):
            if i!=j:
                # temp.append(Distance_polygons(polygons[i],polygons[j],True))
                if Groups[i]==Groups[j]:
                
                    a=Distance_polygons(polygons[i],polygons[j],False)
                    dis_same.append(a)
                else:
                    a=Distance_polygons(polygons[i],polygons[j],True)
                    dis_diff.append(a)
        if len(dis_same)>0:
            sa+=max(dis_same)
        if len(dis_diff)>0:
            di+=min(dis_diff)
    return di-sa



        
        
def overlap_area(A,B): 

    '''trying to caculate the overlap area of two polygons 
    this function works, but it's quite slow
    so it was not used in the whole programm'''

    Ax=A[0::2]
    Ay=A[1::2]
    Bx=B[0::2]
    By=B[1::2]

    temp1=[]
    for i in range(0,len(Ax)):
        temp1.append((Ax[i],Ay[i]))
    temp2=[]
    for i in range(0,len(Bx)):
        temp2.append((Bx[i],By[i]))
    p1=Polygon(*temp1)
    p2=Polygon(*temp2)

    overlap=p1.intersection(p2)

    if len(overlap)==0:
         return 1
    elif len(overlap)==1 and isinstance(overlap[0],Point):
        return 1
    else:
        temp=set([])
        for i in overlap:
            if isinstance(i,Segment):
                temp.add(i.p1)
                temp.add(i.p2)
            elif isinstance(i,Point):
                temp.add(i)
        p1_points=p1.vertices
        p2_points=p2.vertices
        for i in p1_points:
            if p2.encloses_point(i):
                temp.add(i)
        for i in p2_points:
            if p1.encloses_point(i):
                temp.add(i)
        list1=list(temp)
        x0=(list1[0][0]+list1[1][0]+list1[2][0])/3.0
        y0=(list1[0][1]+list1[1][1]+list1[2][1])/3.0
        points=sorted(list1,key=lambda x:atan2(x[1]-y0,x[0]-x0))
        overlap_polygon=Polygon(*points)
        return overlap_polygon.area
    


def Distance_polygons (A,B,flag):
    #flag False means in the same class

    Ax=A[0::2]
    Ay=A[1::2]
    Bx=B[0::2]
    By=B[1::2]
    Aedge=[]
    Bedge=[]
    dist=[]
    


    for i in range(0,len(Ax)-1):
        Aedge.append(tuple((Ax[i],Ay[i],Ax[i+1],Ay[i+1])))
    Aedge.append(tuple((Ax[-1],Ay[-1],Ax[0],Ay[0])))
    for i in range(0,len(Bx)-1):
        Bedge.append(tuple((Bx[i],By[i],Bx[i+1],By[i+1])))
    Bedge.append(tuple((Bx[-1],By[-1],Bx[0],By[0])))
    
    for i in range(0,len(Ax)):
        for j in range(0,len(Bedge)):
            dist.append(Distance_Point2LineSeg(Bedge[j],Ax[i],Ay[i]))
    for i in range(0,len(Bx)):
        for j in range(0,len(Aedge)):
            dist.append(Distance_Point2LineSeg(Aedge[j],Bx[i],By[i]))
    
    if flag:
        if isPolygonsOverlap(A,B):
            # return (-1)*(min(dist))*sqrt(overlap_area(A,B))
            return (-1)*(min(dist))*20
        
        else:
            return min(dist)
    else:
        if isPolygonsOverlap(A,B):
            # return min(dist)*sqrt(overlap_area(A,B))
            return min(dist)*20
        
        else:
            return min(dist)




def isPolygonsOverlap(A,B):

    # return True means overlap

    Ax=A[0::2]
    Ay=A[1::2]
    Bx=B[0::2]
    By=B[1::2]
    Aedge=[]
    Bedge=[]
    
    for i in range(0,len(Ax)-1):
        Aedge.append(tuple((Ax[i],Ay[i],Ax[i+1],Ay[i+1])))
    Aedge.append(tuple((Ax[-1],Ay[-1],Ax[0],Ay[0])))
    for i in range(0,len(Bx)-1):
        Bedge.append(tuple((Bx[i],By[i],Bx[i+1],By[i+1])))
    Bedge.append(tuple((Bx[-1],By[-1],Bx[0],By[0])))

    for i in range(0,len(Aedge)):
        for j in range(0,len(Bedge)):
            if LineIntersection.IfLineIntersect(Aedge[i],Bedge[j]):
                return True
    if(max(Ax)>max(Bx) and min(Ax)<min(Bx) and max(Ay)>max(By) and min(Ay)<min(By)):
        return True
    if(max(Bx)>max(Ax) and min(Bx)<min(Ax) and max(By)>max(Ay) and min(By)<min(Ay)):
        return True
    return False



def Distance_Point2LineSeg (edge, x3,y3): # x3,y3 is the point

    ''' the distance between a point and a line segment'''
    x1=edge[0]
    y1=edge[1]
    x2=edge[2]
    y2=edge[3]
    px = x2-x1
    py = y2-y1
    something = px*px + py*py
    u =  ((x3 - x1) * px + (y3 - y1) * py) / float(something)
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    x = x1 + u * px
    y = y1 + u * py
    dx = x - x3
    dy = y - y3

    dist = math.sqrt(dx*dx + dy*dy)
    return dist


