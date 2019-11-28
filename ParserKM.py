import xml.dom.minidom
import re
from zipfile import ZipFile
import os

class ParserKM():
    """ Класс который пердоставляет методы для получения
    имён и координат из файлов kmz и kml """
    def __init__(self,p):
        filename, file_extension = os.path.splitext(p)
        if file_extension==".kmz":
            kmz = ZipFile(p, 'r')
            kml = kmz.open(kmz.namelist()[0], 'r')
            doc = xml.dom.minidom.parse(kml)
        elif file_extension==".kml":
            doc = xml.dom.minidom.parse(p)
        else:
            raise Exception("Error extension file")
        doc.normalize()
        Xml = doc.documentElement 
        self.Placemark = Xml.getElementsByTagName("Placemark")
    
    def Name(self,s):
        name = []
        for i in range(len(self.Placemark)):
            tags = self.Placemark[i].getElementsByTagName(s)
            if len(tags) != 0 :
                tags = self.Placemark[i].getElementsByTagName("name")[0]
                name.append([self.getText(tags.childNodes),i])
        return name

    def List(self,s,n=None):
        List = []
        if n==None:
            for i in range(len(self.Placemark)):
                tags = self.Placemark[i].getElementsByTagName(s)
                if len(tags) != 0 :
                    tags = self.Placemark[i].getElementsByTagName("coordinates")[0]
                    List.append(self.ParseCord(self.getText(tags.childNodes)))
            return List
        else:
            tags = self.Placemark[n].getElementsByTagName(s)
            if len(tags) != 0 :
                tags = self.Placemark[n].getElementsByTagName("coordinates")[0]
                return self.ParseCord(self.getText(tags.childNodes))

    def PointName(self):
        return self.Name("Point")
    def PointList(self,n=None):
        return self.List("Point",n)

    def LineName(self):
        return self.Name("LineString")
    def LineList(self,n=None):
        return self.List("LineString",n)

    def ParseCord(self,coordinates):
        Pointstr =  re.split("\n\t\t\t\t|\n\t\t\t|,0|,",coordinates)
        Point = []
        k = 0
        for i in range(len(Pointstr)):
            try:
                if k==0: a=[]
                a.append(float(Pointstr[i]))
                k+=1
                if k==2: 
                    Point.append(tuple([a[1],a[0]]))
                    k=0
            except Exception:
                k=0
        if len(Point)==1:
            return Point[0]
        return Point

    def getText(self, nodelist):
        rc = ""
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc = rc + node.data
        return rc



#print(p.PointName())
#print(p.PointList(4))
#print(p.LineName())
#print(p.LineList(7))
#'ВЛ 330 кВ Лукомльская ГРЭС- ПС 330 кВ Борисов'
""" for i in p.LineName():
    print(i) """