import numpy as np
import matplotlib.pyplot as plt
import copy
from ParserKM import ParserKM
import os
import matplotlib.pyplot as plt
plt.figure(figsize=(7,7))


def ApprocsiLine(sp):
    """ Апроксимация методом наименьших квадратов """
    s_xy=0
    s_x=0
    s_y=0
    s_x2=0
    n = len(sp)
    if n<2:return None
    for i in sp:
        s_xy+=i[0]*i[1]
        s_x+=i[0]
        s_y+=i[1]
        s_x2+=i[0]**2
    a = (n*s_xy-s_x*s_y)/(n*s_x2-s_x**2)
    b = (s_y-a*s_x)/n
    return [a,b]
    
def XYcord(kord):
    #print(kord)
    """ Переход от Геодезических координат 
    к плоским прямоугольным координатам """
    B = kord[0]/180*np.pi
    L = kord[1]
    n=int((6+L)/6)
    l=(L-(3+6*(n-1)))/57.29577951
    l2=l**2
    s2=(np.sin(B))**2
    s4=(np.sin(B))**4
    s6=(np.sin(B))**6


    x = 6367558.4968*B-np.sin(2*B)*(16002.8900+66.9607*s2+0.3515*s4-\
        l2*(1594561.25+5336.535*s2+26.790*s4+0.149*s6+\
        l2*(672483.4-811219.9*s2+5420.0*s4-10.6*s6+\
        l2*(278194-830174*s2+572434*s4-16010*s6+\
        l2*(109500-574700*s2+863700*s4-398600*s6)))))
    
    y = (5+10*n)*10**5+l*np.cos(B)*(6378245+21346.1415*s2+107.1590*s4+\
        0.5977*s6+l2*(1070204.16-2136826.66*s2+17.98*s4-11.99*s6+\
        l2*(270806-1523417*s2+1327645*s4-21701*s6+\
        l2*(79690-866190*s2+1730360*s4-945460*s6))))
    return (x,y)



def greatCircleDistance(kord1, kord2,R_eath=6371):
    """ Расчёт растояния по географическим координатам """
    (lat1,lon1) = kord1
    (lat2,lon2) = kord2
    def haversin(x):
        return np.sin((x/180*np.pi)/2)**2 
    return R_eath * (2 * np.arcsin(np.sqrt(haversin(lat2-lat1) 
          +np.cos(lat1/180*np.pi) * np.cos(lat2/180*np.pi) * haversin(lon2-lon1))))

#print(greatCircleDistance((52.0,25.0),(52.5,25.0)))
#print(greatCircleDistance((52.0,25.0),(52.0,25.5)))

def Distance(line):
    """ Расчёт растояния между двумя точками """
    #print(np.sqrt((line[1][0]-line[0][0])**2+(line[1][1]-line[0][1])**2))
    return np.sqrt((line[1][0]-line[0][0])**2+(line[1][1]-line[0][1])**2)

def Perpendicular(line1,line2,diapason=(89,91),trig=False):
    """ Определение перпендикулярности двух прямых """
    [(x1,y1),(x2,y2)] = line1
    [(x3,y3),(x4,y4)] = line2
    if x1==x2:
        u1=np.pi/2*(y2-y1)/abs((y2-y1))
    else:
        u1=np.arctan((y2-y1)/(x2-x1))
    if x3==x4:
        u2=np.pi/2*(y4-y3)/abs((y4-y3))
    else:
        u2=np.arctan((y4-y3)/(x4-x3))
    u=abs(u1-u2)*180/np.pi

    if not trig:
        return diapason[1]>u>diapason[0]
    else:
        return not diapason[1]>u>diapason[0]


def Lies(point,line):
    """ Проверка принадлежности точки отрезку """
    (x4,y4)=point
    [(x1,y1),(x2,y2)]=line
    e1=0.01

    if (abs(x4 - x1)<e1 and abs(y4 - y1)<e1) or (abs(x4 - x2)<e1 and abs(y4 - y2)<e1):
        return True
    else:
        if x1==x2==x4 and (y1<=y4<=y2 or y1>=y4>=y2):
            return True
        elif y1==y2==y4 and (x1<=x4<=x2 or x1>=x4>=x2):
            return True
        elif ((x1<=x4<=x2 or x1>=x4>=x2) and (y1<=y4<=y2 or y1>=y4>=y2)):
            return (x2-x1)/(x4-x1)-(y2-y1)/(y4-y1)<0.01
        else: return False
    


def Intersection(point,line,trig=False):
    """ функция опускания перпендикуляра с точки на линию """
    (x,y)=point
    [(x1,y1),(x2,y2)]=line
    x4=((x2-x1)*(y2-y1)*(y-y1)+x1*pow(y2-y1, 2)+x*pow(x2-x1, 2))/(pow(y2-y1, 2)+pow(x2-x1, 2))
    y4=(y2-y1)*(x4-x1)/(x2-x1)+y1

    
    t1 = Lies((x4,y4),line)
    t3 = Distance([(x,y),(x4,y4)]) <= 300


    if not trig:
        if t1  and t3:
            return (x4,y4)
    else:
        if Distance([(x,y),(x4,y4)]) <= 600:
            return (x4,y4)

def Location(point,line):
    """ Проверка расположения точки относительно линии """
    (x3,y3)=point
    [(x1,y1),(x2,y2)]=line
    D = (x3 - x1) * (y2 - y1) - (y3 - y1) * (x2 - x1)
    if D == 0:
        return 0
    elif D < 0:
        return 1
    elif D > 0:
        return -1

def LiesBadPoint(point,line,h):
    """ Проверка находится ли точка вблизи линии """
    if Lies(point,line):
        return point
    else:
        (x4,y4)=point
        [(x1,y1),(x2,y2)]=line
        if (y4-y1)**2+(x4-x1)**2<=h**2:
            return (x1,y1)
        if (y4-y2)**2+(x4-x2)**2<=h**2:
            return (x2,y2)
        if (x1-h<=x4<=x2+h or x1+h>=x4>=x2-h) and (y1-h<=y4<=y2+h or y1+h>=y4>=y2-h):
            point2 = Intersection(point,line)
            if point2 != None:
                if Distance([point,point2]) < h:
                    return point2

        return False

def VHLine(line):
    """ Получение близких к граничным координат на линии """
    [(x1,y1),(x2,y2)] = line
    d=Distance(line)
    m=0.1/d
    l1=m/(1-m)
    l2=(1-m)/m
    x3=(x1+x2*l1)/(1+l1)
    y3=(y1+y2*l1)/(1+l1)
    x4=(x1+x2*l2)/(1+l2)
    y4=(y1+y2*l2)/(1+l2)
    
    return [(x3,y3),(x4,y4)]

def Psled(rvl,vvl):
    """ Определение взаимного располежения РВЛ и ВВЛ """
    sl =[]
    Doble = set()
    N = set()
    for i in range(1,len(rvl)):#len(rvl)

        a = []
        b = []
        c = []
        for j in range(1,len(vvl)):#len(vvl)
            if not Perpendicular([rvl[i-1],rvl[i]],[vvl[j-1],vvl[j]],diapason=(70,130),trig=True):
                continue
            
            rvl1 = VHLine([rvl[i-1],rvl[i]])
            Yz = [rvl1[0],rvl1[1],vvl[j-1],vvl[j]]
            Line = [[vvl[j-1],vvl[j]],[vvl[j-1],vvl[j]],rvl1, rvl1]

            Base = [a,b,c,c]
            K = [-1,-1,-1,-1]
            Rasch = []

            for k in range(4):
                if k<2:
                    Rasch.append(Intersection(Yz[k],Line[k],trig=True))
                    if Rasch[k] == None:
                        if Yz[k] not in N:
                            N.add(Yz[k])
                            Doble.add((Yz[k],None))
                            #Base[k].append([Yz[k],None])
                    else:
                        if (Yz[k],Rasch[k]) not in Doble:
                            N.add(Yz[k])                  
                            obr = Intersection(Rasch[k],Line[2],True)

                            if obr != None:
                                cos_alf = Distance([obr,Rasch[k]])/Distance([Rasch[k],Yz[k]])
                                alf=np.arccos(cos_alf)
                                #bet=np.pi/2-alf
                                #sin_bet = np.sin(bet) 
                                if cos_alf != 1:
                                    x3_1 = (Yz[k][0]/np.tan(np.pi/2)+Rasch[k][0]/np.tan(alf)+Yz[k][1]-Rasch[k][1])/(1/np.tan(alf)+1/np.tan(np.pi/2))
                                    y3_1 = (Yz[k][1]/np.tan(np.pi/2)+Rasch[k][1]/np.tan(alf)+Rasch[k][0]-Yz[k][0])/(1/np.tan(alf)+1/np.tan(np.pi/2))

                                    x3_2 = (Yz[k][0]/np.tan(np.pi/2)+Rasch[k][0]/np.tan(alf)+Rasch[k][1]-Yz[k][1])/(1/np.tan(alf)+1/np.tan(np.pi/2))
                                    y3_2 = (Yz[k][1]/np.tan(np.pi/2)+Rasch[k][1]/np.tan(alf)+Yz[k][0]-Rasch[k][0])/(1/np.tan(alf)+1/np.tan(np.pi/2))

                                    obr1 = Intersection((x3_1,y3_1),Line[2])
                                    obr2 = Intersection((x3_2,y3_2),Line[2])
                                    #print(obr1,obr2)
                                        
                                    if obr1 != None and Perpendicular([Yz[k],(x3_1,y3_1)],Line[2]) and Lies((x3_1,y3_1),Line[k]):#
                                        Doble.add((Yz[k],(x3_1,y3_1)))
                                        Base[k].append([Yz[k],(x3_1,y3_1)])
                                        K[k]+=1
                                    if obr2 != None and Perpendicular([Yz[k],(x3_2,y3_2)],Line[2]) and Lies((x3_2,y3_2),Line[k]):#
                                        Doble.add((Yz[k],(x3_2,y3_2)))
                                        Base[k].append([Yz[k],(x3_2,y3_2)])
                                        K[k]+=1

                                else:
                                    if Lies(Rasch[k],Line[k]):
                                        Doble.add((Yz[k],Rasch[k]))
                                        Base[k].append([Yz[k],Rasch[k]])
                                        K[k]+=1

                            else:
                                if Lies(Rasch[k],Line[k]):
                                    Doble.add((Yz[k],Rasch[k]))
                                    Base[k].append([Yz[k],Rasch[k]])
                                    K[k]+=1

                else:
                    Rasch.append(Intersection(Yz[k],Line[k]))
                    if Rasch[k] != None and (Rasch[k],Yz[k]) not in Doble:
                        N.add(Rasch[k])
                        Doble.add((Rasch[k],Yz[k]))
                        Base[k].append([Rasch[k],Yz[k]])
            #print(Rasch)
        
        l = [[Distance([c[j][0],rvl[i-1]]),j] for j in range(len(c))]
        ls = sorted(l)

        sl+=a
        for j in range(len(ls)):
            sl.append(c[ls[j][1]])
        sl+=b

    #print(sl)
    # Переводим координаты в километры от начала линии
    d_v = set()
    k=0
    km=0
    rvl_ych = []
    vvl_ych = []
    a = []
    b = []
    la=0
    #print(len(sl))
    #print("---------")
    for i in range(1,len(rvl)):
        for k in range(len(sl)):
            if Lies(sl[k][0],[rvl[i-1],rvl[i]]) and k not in d_v:
                lr = Location(sl[k][1],[rvl[i-1],rvl[i]])
                h = Distance([sl[k][0],sl[k][1]])
                d = Distance([rvl[i-1],sl[k][0]])
                a.append([d+km,h*lr])
                b.append(sl[k][1])
                d_v.add(k)
   
        km+=Distance([rvl[i-1],rvl[i]])

        if la==len(a) and len(a)!=0:
            rvl_ych.append(a)
            vvl_ych.append(b)
            la=0
            a=[]
            b=[]
        else:
            la=len(a)
        
    if la==len(a) and len(a)!=0:
        rvl_ych.append(a)
        vvl_ych.append(b)
        la=0
        a=[]
        b=[]
    else:
        la=len(a)
    
    
    return rvl_ych, vvl_ych, km, sl

def Stabilization(rvl,vvl):
    """ Обработка и коректировка участков сближения РВЛ и ВВЛ """

    # Убираем точки перевёртыши
    rvl1 = []
    vvl1 = []
    for i in range(len(rvl)):
        k = len(rvl[i])
        if k==2:
            if not ((rvl[i][0][1]>0 and rvl[i][1][1]<0) or (rvl[i][0][1]<0 and rvl[i][1][1]>0)):
                rvl1.append(rvl[i])
                vvl1.append(vvl[i])
        else:
            for j in range(k-2):
                if (rvl[i][j][1]>0 and rvl[i][j+1][1]<0 and rvl[i][j+2][1]>0) or \
                    (rvl[i][j][1]<0 and rvl[i][j+1][1]>0 and rvl[i][j+2][1]<0):
                    rvl[i][j+1][1]*=-1
            rvl1.append(rvl[i])
            vvl1.append(vvl[i])

    # Убираем повторяющиеся точки
    rvl2 = []
    vvl2 = []
    for i in range(len(rvl1)):
        srez = [0]
        k = len(rvl1[i])
        if k>=2:
            if rvl1[i][0][0] == rvl1[i][1][0]:
                del rvl1[i][0]
                del vvl1[i][0]
            k = len(rvl1[i])
            if k>=2:
                if rvl1[i][k-2][0] == rvl1[i][k-1][0]:
                    del rvl1[i][k-1]
                    del vvl1[i][k-1]
                j=1
                while j<len(rvl1[i]) and len(rvl1[i])>=2:
                    if rvl1[i][j-1][0] == rvl1[i][j][0]:
                        if abs(rvl1[i][j-1][1]-rvl1[i][j][1])>10:
                            srez.append(j)
                            j+=1
                        else:
                            srd = (rvl1[i][j-1][1]+rvl1[i][j][1])/2
                            rvl1[i][j-1][1] = srd
                            del rvl1[i][j]
                            del vvl1[i][j]
                    else: 
                        j+=1
                k = len(rvl1[i])
                srez.append(k)
                for j in range(1,len(srez)):
                    if srez[j]-srez[j-1]>=2:
                        rvl2.append(rvl1[i][srez[j-1]:srez[j]])
                        vvl2.append(vvl1[i][srez[j-1]:srez[j]])

    # Деление
    rvl3 = []
    vvl3 = []
    for i in range(len(rvl2)):
        srez = [0]
        try:
            lr_st = abs(rvl2[i][0][1])/rvl2[i][0][1]
        except ZeroDivisionError:
            lr_st = 1
        for j in range(1,len(rvl2[i])):
            try:
                lr_n = abs(rvl2[i][j][1])/rvl2[i][j][1]
            except ZeroDivisionError:
                lr_n = lr_st
            
            if (lr_st>0 and lr_n<0) or (lr_st<0 and lr_n>0):
                srez.append(j)
            lr_st=lr_n
        
        srez.append(len(rvl2[i]))

        for j in range(1,len(srez)):
            if srez[j]-srez[j-1]>=2:
                rvl3.append(rvl2[i][srez[j-1]:srez[j]])
                vvl3.append(vvl2[i][srez[j-1]:srez[j]])


    # Разбиение на апроксимируемые участки
    rvl4 = []
    vvl4 = []
    for j in range(len(rvl3)):
        v = []
        ind=[]
        k = -1
        for i in range(1,len(rvl3[j])):
            if rvl3[j][i][0]==rvl3[j][i-1][0]:
                continue
            #[a,b]=ApprocsiLine([rvl3[j][i-1],rvl3[j][i]])
            # Замеряем градус отклонения от оси Х
            b= 1 if rvl3[j][i][1] >= rvl3[j][i-1][1] else -1
            a = np.arccos((rvl3[j][i][0]-rvl3[j][i-1][0])/Distance([rvl3[j][i-1],rvl3[j][i]]))*b
            v.append(a)
            ind.append(i)
            k+=1
        if rvl3[j][i][0]==rvl3[j][i-1][0]: 
            ind[k]=i
        k=0
        kn=0
        srez = []
        while k<len(v):
            for i in range(k,len(v)):
                a_sr = abs(sum(v[k:i+1])/(i+1-k))
                a_m = max(abs(max(v[k:i+1])),abs(min(v[k:i+1])))
                if abs(a_sr-a_m)<=15/180*np.pi: # граница отклонения 5 градусов
                    1
                else:
                    srez.append([kn,ind[i]])
                    k=i
                    kn=k
                    break
                if i == len(v)-1:
                    k=len(v)
                    srez.append([kn,ind[k-1]+1])

        for i in srez:
            if i[1]-i[0]>=2:
                rvl4.append(rvl3[j][i[0]:i[1]])
                vvl4.append(vvl3[j][i[0]:i[1]])
    
    rvl6 = copy.deepcopy(rvl4)
    vvl6 = copy.deepcopy(vvl4)
    # Апроксимация обособленых участков
    rvl5 = []
    vvl5 = []
    for i in range(len(rvl4)):
        k = len(rvl4[i])
        if k>2:
            [a,b]=ApprocsiLine(rvl4[i])
            p1 = [rvl4[i][0][0],a*rvl4[i][0][0]+b]
            p2 = [rvl4[i][k-1][0],a*rvl4[i][k-1][0]+b]
            rvl5.append([p1,p2])
            vvl5.append([vvl4[i][0],vvl4[i][k-1]])
        else:
            rvl5.append(rvl4[i])
            vvl5.append(vvl4[i])

    return rvl5, rvl6, vvl5, vvl6

def KnOp(vl,Op_point,h):
    """ Перевод координат опор в км от начала линии """
    km = 0
    Op = []
    d = set()
    for i in range(1,len(vl)):
        for j in range(len(Op_point)):
            if j not in d:
                point2 = LiesBadPoint(Op_point[j][1],[vl[i-1],vl[i]],h)
                if point2:
                    Op.append([Op_point[j][0],km+Distance([vl[i-1],point2])])
                    d.add(j)
            else:
                continue

        km+=Distance([vl[i-1],vl[i]])
    
    return Op

def KmVVL(sp1,vvl):
    """ Перевод координат пересечения ВВЛ с РВЛ в км от начала линии """
    sp = copy.deepcopy(sp1)
    km=0
    d=set()
    for i in range(len(vvl)):
        for j in range(len(sp)):
            for k in range(len(sp[j])):
                if (j,k) not in d:
                    if Lies(sp[j][k],[vvl[i-1],vvl[i]]):
                        sp[j][k]=km+Distance([vvl[i-1],sp[j][k]])
                        d.add((j,k))

        km+=Distance([vvl[i-1],vvl[i]])
    return sp

def NumOpor(sp,Op):
    """ Перевод  значений км в значения номеров опор"""
    #sp1 = copy.deepcopy(sp)
    sp1 = []
    for i in range(len(sp)):
        a=[]
        for j in range(len(sp[i])):
            a.append(copy.deepcopy(sp[i][j]))
        sp1.append(a)

    d = set()
    for i in range(len(sp1)):
        for j in range(len(sp1[i])):
            for k in range(1,len(Op)):
                s = str(type(sp1[i][j]))
                if (i,j) not in d:
                    if s == "<class 'int'>" or s == "<class 'float'>" or s=="<class 'numpy.float64'>":
                        t1 =  1 if Op[k-1][1]<=sp1[i][j]<=Op[k][1] else 0
                        t2 = -1 if Op[k-1][1]>=sp1[i][j]>=Op[k][1] else 0
                        t3 = 1 if Op[k-1][0]<Op[k][0] else -1
                        if t1==1 or t2==-1:
                            Pr = abs(Op[k][1]-Op[k-1][1])/abs(Op[k][0]-Op[k-1][0])
                            O = int(round(abs(sp1[i][j]*(t1+t2)-Op[k-1][1]*(t1+t2))/Pr*t3+Op[k-1][0]))
                            sp1[i][j] = O
                            d.add((i,j))

                    elif s == "<class 'list'>":
                        t1 =  1 if Op[k-1][1]<=sp1[i][j][0]<=Op[k][1] else 0
                        t2 = -1 if Op[k-1][1]>=sp1[i][j][0]>=Op[k][1] else 0
                        t3 = 1 if Op[k-1][0]<Op[k][0] else -1
                        if t1==1 or t2==-1:
                            Pr = abs(Op[k][1]-Op[k-1][1])/abs(Op[k][0]-Op[k-1][0])
                            O = int(round(abs(sp1[i][j][0]*(t1+t2)-Op[k-1][1]*(t1+t2))/Pr*t3+Op[k-1][0]))
                            sp1[i][j][0] = O
                            d.add((i,j))
    
    return sp1

def CheckVL(name,sp):
    """ Проверяем наличие искомой ВЛ """
    j = None
    for i in range(len(sp_vl)):
        if name in sp_vl[i][0]:
            j=i
            break
    return j

def Otp(rvl, sp_otp,sp_vl,BD_VL):
    """ Определяем места премыкания отпаек к расматриваемой ВЛ """
    l = []
    n = []
    nch_otp =[]
    for i in range(len(sp_otp)):
        a=[]
        j = CheckVL(sp_otp[i],sp_vl)
        if j != None:
            o_geo = BD_VL.LineList(sp_vl[j][1])
            o1 = XYcord(o_geo[0]) 
            o2 = XYcord(o_geo[len(o_geo)-1])
            km = 0
            for k in range(1,len(rvl)):
                t1 = Intersection(o1,[rvl[k-1],rvl[k]])
                t2 = Intersection(o2,[rvl[k-1],rvl[k]])
                if t1 != None and t2 !=None:
                    if Distance([t1,o1])<Distance([t2,o2]):
                        t=t1
                        o=o1
                        nch=1
                    else:
                        t=t2
                        o=o2
                        nch=-1
                    
                elif t1 != None:
                    t=t1
                    o=o1
                    nch=1
                elif t2 !=None:
                    t=t2
                    o=o2
                    nch=-1
                else: t=None

                if t != None:
                    if Distance([t,o])<10:
                        dst = km+Distance([rvl[k-1],t])
                        a.append(int(round(dst)))
                else:
                    dd1=Distance([rvl[k-1],o1])
                    dd2=Distance([rvl[k],o1])
                    dd3=Distance([rvl[k-1],o2])
                    dd4=Distance([rvl[k],o2])
                    if dd1<10 or dd3<10:
                        nch = 1
                        dst = km
                        a.append(int(round(dst)))
                    elif dd2<10 or dd4<10:
                        nch = -1
                        dst = km+Distance([rvl[k-1],rvl[k]])
                        a.append(int(round(dst)))
                    
                
                km+=Distance([rvl[k-1],rvl[k]])

            if a != []:
                l.append(int(round(sum(a)/len(a))))
                n.append(sp_otp[i])
                nch_otp.append(nch)
    
    return l, n, nch_otp

def Plots(rvl,vvl,sl):
    """ Построение графика сближения ВЛ на месности """
    x1 = [i[0] for i in rvl]
    y1 = [i[1] for i in rvl]

    x2 = [i[0] for i in vvl]
    y2 = [i[1] for i in vvl]

    Deistv = plt.figure(figsize=(7,7))
    ax=Deistv.add_subplot(111)

    ax.plot(x1, y1,"r",label="rvl")
    ax.plot(x2, y2,"g",label="vvl")

    for i in sl:
        ax.plot((i[0][0],i[1][0]), (i[0][1],i[1][1]),"y")

    ax.grid(True) # Включаем сетку
    ax.set_xlabel('X') # Подпись оси х
    ax.set_ylabel('Y') # Подпись оси у
    ax.legend(frameon=False) # Выводим легенду графика

    mx_x = -float('inf')
    mn_x = float('inf')

    mx_y = -float('inf')
    mn_y = float('inf')

    for i in rvl+vvl:
        if i[0]>mx_x: mx_x= i[0]
        if i[1]>mx_y: mx_y= i[1]
        if i[0]<mn_x: mn_x= i[0]
        if i[1]<mn_y: mn_y= i[1]

    m = max(mx_x-mn_x,mx_y-mn_y)
    ax.axis([mn_x,mn_x+m,mn_y,mn_y+m])

    Deistv.show()

def OpenFile():
    """ Чтения зандания на поиск из файла txt """
    #fname = str(input())
    fname = "inp.txt"
    m =[]
    sp_otp = []
    with open(fname, 'r') as f:
        for line in f:
            m.append(str(line))
    i=0
    while i<len(m):
        if m[i] == "Файл карт:\n":
            cart = m[i+1].replace("\n","").strip()
            i+=2
        elif m[i] == "Название линии:\n":
            name_napr, name = m[i+1].split(":")
            name = name.replace("\n","").strip()
            name_napr = int(name_napr)
            i+=2
        elif m[i] == "Название отпайки:\n":
            i+=1
            continue
        else:
            sp_otp.append(m[i].replace("\n","").strip())
            i+=1

    return cart, name, name_napr, sp_otp
            


cart, name, name_napr, sp_otp= OpenFile()    

            
BD_VL = ParserKM("google_maps_kmz/"+cart)



#name = "ВЛ 35 кВ ПС 110 кВ Толочин - ПС 35 кВ Литвяки"
#sp_otp = ["ВЛ 35 кВ отпайка на ПС 35 кВ Толочин - 2"]

i = CheckVL(name,sp_vl)
print(i)

print(sp_vl[i][0])
txt = open("Results/"+sp_vl[i][0]+'.txt','w')

print("----------")
txt.write(sp_vl[i][0]+'\n')
txt.write("----------"+'\n')

rvl_geo = BD_VL.LineList(sp_vl[i][1])#[::-1]

rvl = [XYcord(i) for i in rvl_geo]
if name_napr == -1:
    rvl.reverse()


name_vvl = []
vl_ych = []



otp, name_otp, nch_otp = Otp(rvl, sp_otp,sp_vl,BD_VL)

for j in range(0,len(sp_vl)):#len(sp_vl)
    vvl_geo = BD_VL.LineList(sp_vl[j][1])#sp_vl[j][1]353
    
    try:
        vvl = [XYcord(i) for i in vvl_geo]
    except Exception:
        continue

    if sp_vl[j][1] == sp_vl[i][1]:
        continue

    
    
    rvl_ych, vvl_ych, km, sl = Psled(rvl,vvl)


    if rvl_ych != []:   
        r5,r6,v5,v6 = Stabilization(rvl_ych,vvl_ych)
        if r6 !=[]:
            name_vvl.append(sp_vl[j][0])
            vl_ych.append(r5)



d = {}

for i in range(len(vl_ych)):
    for j in range(len(vl_ych[i])):
        n=int(round(vl_ych[i][j][0][0]))
        k=int(round(vl_ych[i][j][1][0]))
        if n not in d:
            d[n]=0
        if k not in d:
            d[k]=0


d[0]=0
d[int(round(km))]=0

for i in range(len(otp)):
    d[otp[i]]=0

s_otp = set(otp)

l = sorted(list(d.keys()))

m=0
n=m+1
while m<len(l):
    while n<len(l):
        if l[n]-l[m]<float(sp_otp[3]) and n!=len(l)-1 and l[n] not in s_otp:
            del l[n]
        elif l[n]-l[m]<float(sp_otp[3]) and n==len(l)-1 and l[m] not in s_otp:
            del l[m]
        elif l[n]-l[m]<float(sp_otp[3]) and l[n] in s_otp and l[m] not in s_otp:
            del l[m]   
        else:
            break
    m+=1
    n=m+1


for i in range(len(l)):
    d[l[i]] = i+1

d2 = {}
l2 = sorted(list(d.keys()))


for i in range(len(l2)):
    k=float("inf")
    ind=0
    for j in range(len(l)):
        dl = abs(l[j]-l2[i])
        if dl<k:
            #print("test")
            k=dl
            ind = j
    d2[l2[i]] = d[l[ind]]



n_y_otp=[]
for i in range(1,len(l)):
    print(d[l[i-1]],d[l[i]],l[i]-l[i-1])
    txt.write(str(d[l[i-1]])+' '+str(d[l[i]])+' '+str(l[i]-l[i-1])+'\n')
    if l[i-1] in s_otp:
        n_y_otp.append([d[l[i-1]],name_otp[otp.index(l[i-1])],nch_otp[otp.index(l[i-1])]])
if l[i] in s_otp:
    n_y_otp.append([d[l[i]],name_otp[otp.index(l[i])],nch_otp[otp.index(l[i])]])

print("-----------")

for i in n_y_otp:
    print(i)
    txt.write(str(i)+'\n')

print("-----------")


vl_ych2 = []
name_vvl2 = []
g = -1
for i in range(len(vl_ych)):
    g+=1
    vl_ych2.append([])
    h=-1
    for j in range(len(vl_ych[i])):
        n=d2[int(round(vl_ych[i][j][0][0]))]
        k=d2[int(round(vl_ych[i][j][1][0]))]
        if k-n != 0:
            h+=1 
            vl_ych2[g].append(np.array(vl_ych[i][j]))
            vl_ych2[g][h][0][0]=n
            vl_ych2[g][h][1][0]=k

    
    if vl_ych2[g]==[]:
        del vl_ych2[g]
        g-=1
    else:
        name_vvl2.append(name_vvl[i])


for i in range(len(name_vvl2)):
    print(name_vvl2[i])
    txt.write(name_vvl2[i]+'\n')
    for j in range(np.shape(vl_ych2[i])[0]):
        print(vl_ych2[i][j])
        print("-------")
        txt.write(str(vl_ych2[i][j])+'\n')
        txt.write("----------"+'\n')

txt.close()
