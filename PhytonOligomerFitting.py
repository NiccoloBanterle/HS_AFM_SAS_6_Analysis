#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:58:36 2018

@author: banterle
"""
from joblib import Parallel, delayed
from pathlib import Path
import multiprocessing
import timeit
import cv2
import numpy as np
import pprint
import matplotlib.pyplot as plt
import os
from scipy.spatial import distance
import math
import scipy.optimize as sciOpt
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
import json
import pickle
import re
from collections import Counter
from PIL import Image
import h5py   
import pandas as pd
from natsort import natsorted
import operator
class SAS6:
    
    def __init__(self, nativeAngle = 0.666, nativeLength = 11, tempLineWidth = 2, ringclosed = False):
        
        self.elementAngles = []
        self.elementVectors = []
        self.elementPos = []      
        self.templateLineWidth = tempLineWidth
        self.templateDilation = 12
        self.nativeAngle = nativeAngle
        self.nativeLength = nativeLength
        self.closed = ringclosed
        self.matchLoc = (0,0)
        self.classnum = []
    def SetSavePath(self,savePath):
        self.savePath = savePath
    def loadImageData(self, imagePath):
        self.imageData = cv2.imread(imagePath, cv2.IMREAD_GRAYSCALE)
        self.imageSize = self.imageData.shape
    def loadSegmentation(self, fileName=False):
        self.segmentation = ImageContainer(fileName)
        self.classnum = self.segmentation.getClasses()
        
    def loadDirectImageData(self, image):    
        self.imageData = image
        self.imageSize = self.imageData.shape
    def showImageData(self):
        try:
            cv2.imshow('Stored image',self.imageData)
            cv2.waitKey(0)
        except AttributeError:
            print('Image data not loaded.')
        else:
            pass
        
               
    def showTemplate(self):
        try:
            cv2.imshow('Template image',self.template)
            cv2.waitKey(0)
        except AttributeError as e:
            print('Image data not loaded: {0}'.format(e))
        else:
            pass
    
    def showTemplateOverlay(self):
        try:
            rows,cols= self.template.shape
            overlay=cv2.addWeighted(self.imageData[self.matchLoc[0]:self.matchLoc[0]+rows, self.matchLoc[1]:self.matchLoc[1]+cols],0.5,self.template,0.5,0)
            totalOverlay = self.imageData.copy()
            totalOverlay[self.matchLoc[0]:self.matchLoc[0]+rows, self.matchLoc[1]:self.matchLoc[1]+cols] = overlay
            cv2.imshow('Template overlay',totalOverlay)
            cv2.waitKey(0)
        except AttributeError as e:
            print('Image data not loaded: {0}'.format(e))
        else:
            pass
            
    def makeTemplateFromPoints(self):
        drawOffset = np.min(self.elementPos,axis=0)
        templateImageSize = np.max(self.elementPos,axis=0)-drawOffset+np.ones(2)*(self.templateDilation+self.templateLineWidth)
        self.template = np.zeros((int(templateImageSize[1]), int(templateImageSize[0]), 1), np.uint8)
        drawOffset = drawOffset - np.ones(2)*(self.templateDilation+self.templateLineWidth)/2
        
        pointNo = len(self.elementPos)
        
        if self.closed:
            pointInds = range(pointNo)
        else:
            pointInds = range(pointNo-1)
        
        for pointInd in pointInds:
            
            #print(tuple(self.elementPos[pointInd]))
            self.template = cv2.line(self.template, tuple([int(x) for x in self.elementPos[pointInd]-drawOffset]), tuple([int(x) for x in self.elementPos[(pointInd+1) % pointNo]-drawOffset]), (255), self.templateLineWidth)
            
                
        if(self.template.shape[0]>(self.imageSize[0]-1) or self.template.shape[1]>(self.imageSize[1]-1)):
                self.template = self.template[ 0: min(self.template.shape[0], self.imageSize[0]-1),0: min(self.template.shape[1], self.imageSize[1]-1)]
        self.template = cv2.GaussianBlur(self.template, (11,11), 0)
        
        #plt.imshow(self.template)
    def matchTemplate(self):
        correl = cv2.matchTemplate(self.imageData, self.template, cv2.TM_CCORR_NORMED)
        min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(correl)
        if max_val==0:
            max_val=0.00001
        return (1/max_val), max_loc #cv2.sumElems(self.template)[0]
        
    def setNativeAngle(self, nativeAngle):
        self.nativeAngle = nativeAngle
     
    def setElementAnglesFromDeltas(self, angleDeltas, initialAngle = 0):
        # define set of angles for ring-like element
        self.elementAngles = np.cumsum([self.nativeAngle+x for x in angleDeltas]) + initialAngle
    def setPointsFromAngles(self):
        # calculate point vectors as function of angles.
        self.elementVectors = [[self.nativeLength*np.cos(x), self.nativeLength*np.sin(x)] for x in self.elementAngles]
        # integrate up to point coordinates
        self.elementPos = np.cumsum(self.elementVectors,axis=0)
        self.elementPos = np.insert(self.elementPos, 0, [0,0], axis=0)
        
    def makeElementsFromDeltas(self, angleDeltas, initialAngle = 0):
        if self.closed: 
            self.setNativeAngle(2*math.pi/(len(angleDeltas)+1))
            self.setElementAnglesFromDeltas(np.zeros(len(angleDeltas)), initialAngle)
            self.setPointsFromAngles() #does it make sense?
            for i in range(3,len(angleDeltas)):
                if not self.perturbElementAngle(i,angleDeltas[i]):
                    return False
            #self.setPointsFromAngles()    
            ## for angle in angledelta call perturbelements
        else:
            self.setNativeAngle(2*math.pi/9)
            self.setElementAnglesFromDeltas(angleDeltas, initialAngle)
            self.setPointsFromAngles()
        return True
        
    def ccw(self,A,B,C):
        return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])

    def intersection(self,seg1, seg2):
        return self.ccw(seg1[0],seg2[0],seg2[1]) != self.ccw(seg1[1],seg2[0],seg2[1]) and self.ccw(seg1[0],seg1[1],seg2[0]) != self.ccw(seg1[0],seg1[1],seg2[1])
        
    def checkPointIntersections(self):
        numPts = len(self.elementPos)
        if numPts > 3:
            # check intersection of current line with all higher ones
            for pointInd in range(numPts-4):
                for pointWalk in range(pointInd+2,numPts-1):
                    if self.intersection([self.elementPos[x] for x in [pointInd, pointInd+1]], [self.elementPos[x] for x in [pointWalk, pointWalk+1]]):
                        return True
            return False
        
        
    def showAngles(self):
        plt.scatter(*zip(*self.elementVectors))
        plt.show()
        
    def showPoints(self):
        plt.plot(*zip(*self.elementPos))
        plt.show()
        
    def getAngleBetweenPoints(self,ptInd1,ptInd2):
        return math.atan2(self.elementPos[ptInd2][1]-self.elementPos[ptInd1][1], self.elementPos[ptInd2][0]-self.elementPos[ptInd1][0])
    
    def getLengthBetweenPoints(self,ptInd1,ptInd2):
        return ( (self.elementPos[ptInd2][1]-self.elementPos[ptInd1][1])**2 + (self.elementPos[ptInd2][0]-self.elementPos[ptInd1][0])**2 )**0.5
    
    def perturbElementAngle(self, element, angleDeltaNew):
        numPts = len(self.elementPos)
        
        elements = [(element-x)% (numPts) for x in range(4)]
        
        #print(elements)
        
        rotationAngle = self.getAngleBetweenPoints(elements[0], elements[3])
        #print(rotationAngle*180/math.pi)
        newAngles = np.zeros(4)
        newAngles[0] = self.getAngleBetweenPoints(elements[0], elements[1])-rotationAngle+angleDeltaNew

        #print([180/math.pi*x for x in newAngles])
        
        C=1/(2*self.nativeLength**2);
        D=C*(2*self.nativeLength**2);
        F=1;

        #print('C: %g D: %g F: %g' % (C, D, F) )


        r =  self.getLengthBetweenPoints(elements[0], elements[3]) - self.nativeLength*math.cos(newAngles[0]);
        s = self.nativeLength*math.sin(newAngles[0]);

        #print('r: %g s: %g  distance: %g' % (r,s, self.getLengthBetweenPoints(elements[0], elements[3])) )

        E=D-C*(r**2+s**2);

        try: 
            delta = math.acos(E);
        except ValueError:
            return False # if problem is infeasible, quit and return an error

        #print('Delta: %f' % (delta*180/math.pi) )

        A=(F-E)/math.sin(delta);
        B=s/r;

        newAngles[1] = math.atan2(1-A*B,A+B);
        newAngles[3] = newAngles[1] + delta;
        newAngles[2] = newAngles[3] - math.pi;
        
        #print([180/math.pi*x for x in self.elementAngles[elements]])
        #print([180/math.pi*(x) for x in newAngles])
        
        self.elementAngles[elements[3]] = newAngles[3] + rotationAngle
        self.elementAngles[elements[2]] = math.pi + newAngles[1] + rotationAngle
        self.elementAngles[elements[1]] = math.pi + newAngles[0] + rotationAngle
        
        self.setPointsFromAngles()
        
        return True
        
    def runOptimization(self, x0 = np.zeros(5)):
        minimizer_kwargs = {"method":"BFGS", "args": ()}
        mybounds = MyBounds()
        mytakestep = MyTakeStep()
        ret = sciOpt.basinhopping(self.costFun, x0, minimizer_kwargs=minimizer_kwargs, niter=20000, T=8, stepsize=0.3, interval = 1000, niter_success=5000 , accept_test=mybounds, take_step=mytakestep) #T=0.8
        return ret
        #myRingFit.overlay_image_mask()
    #def runOptimizationClosed(self, x0 = np.zeros(5)):
        #minimizer_kwargs = {"method":"BFGS", "args": ()}
        #mybounds = MyBounds()
        #ret = sciOpt.basinhopping(self.costFunClosed, x0, minimizer_kwargs=minimizer_kwargs, niter=5000, T=1, stepsize=.5, interval = 15, niter_success=200 , accept_test=mybounds)
        #return ret        
    def costFun(self, x, *args):
        if not self.makeElementsFromDeltas(x[1:], initialAngle = x[0]):
            return 100000
        self.makeTemplateFromPoints()
        #plt.imshow(self.template)
        optiVal, max_loc = self.matchTemplate()
        self.matchLoc = max_loc
        if (self.closed and 1/optiVal>0.97): 
            optiVal=1/(1+optiVal*0.02)
        return optiVal
    #def costFunClosed(self, x, *args):
        #self.makeElementsFromDeltas(x[1:], initialAngle = x[0])
        #self.
        #self.makeTemplateFromPoints()
        #optiVal, max_loc = self.matchTemplate()
        #self.matchLoc = max_loc
        #return optiVal    
    def overlay_image_mask(self):
        img=self.imageData
        img_overlay=self.template 
        #img=np.asarray(Image.fromarray(img).convert('RGBA'))
        img.setflags(write=1)
        x = self.matchLoc[0]
        y = self.matchLoc[1]
        # Image ranges
        x1, x2 = max(0, x), min(img.shape[1], x + img_overlay.shape[1])
        y1, y2 = max(0, y), min(img.shape[0], y + img_overlay.shape[0])

        # Overlay ranges
        x1o, x2o = max(0, -x), min(img_overlay.shape[1], img.shape[1] - x)
        y1o, y2o = max(0, -y), min(img_overlay.shape[0], img.shape[0] - y)
        # Exit if nothing to do
        if y1 >= y2 or x1 >= x2 or y1o >= y2o or x1o >= x2o:
            return
        img[:][:]=0
        img[y1:y2, x1:x2] = img_overlay[y1o:y2o, x1o:x2o]
        #alpha = img_overlay[y1o:y2o, x1o:x2o]/img_overlay.max()
        #alpha_inv = 1.0 - alpha
    
        #img[y1:y2, x1:x2,0] = alpha * img_overlay[y1o:y2o, x1o:x2o]/(img_overlay.max()/255) + alpha_inv * img[y1:y2, x1:x2,0]
        #plt.imshow(img)
        return(img)
    def ImportImageX(self,imageno): #,v_alpha,v_Cmap='Reds'
        
        self.imageData  = self.segmentation.getImage(imageno)
        #print(image.dtype)
        #print(image.shape)
        #plt.imshow(self.imageData)
        #time=(f['table'][imageno][1])
    #    PClass=(f['table'][imageno][4])
        return (imageno,self.imageData)
    def SaveBestFitResult(self,results,i):
        self.closed = (results[0]>9)
        self.makeElementsFromDeltas(results[2][1:], initialAngle = results[2][0])
        self.makeTemplateFromPoints()
        optiVal, max_loc = self.matchTemplate()
        self.matchLoc = max_loc
        imgname = 'img'+str(i)+'.png'
        tmpname = "tmp"+str(i)+'.png'
        ImfileToWriteTo =  self.savePath / imgname
        TmpFileToWriteTo = self.savePath / tmpname
        cv2.imwrite(str(ImfileToWriteTo),self.imageData)
        cv2.imwrite(str(TmpFileToWriteTo),self.overlay_image_mask())

    def processSingleImage(self,i):
            j = self.classnum[i]
            p=j
            current = [-1,0,0]
            while p!=0 or p!=15:
                if(j==p):
                    self.ImportImageX(i)
                    prev = self.BestNFit(self.imageData,j,p,i)
                    if (prev[0]==j):
                        self.SaveBestFitResult(prev,i)
                        print(i)
                        return prev
                    else:
                        if prev[0]<j:
                            if prev[0]>1:
                                p-=1
                            else:
                                self.SaveBestFitResult(prev,i)
                                print(i)
                                return prev
                        else:
                            if prev[0]<14:
                                p+=1
                            else:
                                self.SaveBestFitResult(prev,i)
                                print(i)
                                return prev
                else:
                    #self.ImportImageX(i)
                    current = self.BestNFit(self.imageData,j,p,i)
                    if(current[1]<prev[1]):
                        print(i)
                        self.SaveBestFitResult(prev,i)
                        return prev
                    else:
                        if p>j: 
                            if current[0]<14:
                                prev = current
                                p+=1
                            else:
                                self.SaveBestFitResult(current,i)
                                print(i)
                                return current
                        else:   
                            if current[0]>1:
                                prev = current
                                p-=1
                            else:
                                self.SaveBestFitResult(current,i)
                                print(i)
                                return current          
    def BestNFit(self,image,n,p,imageno):
        nlist=[]
        if n==p:
            if n<=10:
                for j in range(max(int(n-1),1),min(int(n+2),11)): 
                    nlist.append(j)
                    if (j>=7 and j<=10):
                        nlist.append(j+4) ###might be a three           
            if n>10:
                for j in range(max(int(n-1),11),min(int(n+2),15)):
                    nlist.append(j)
                    nlist.append(j-4)###might be a three
        else:
            if p<n:
                    nlist.append(p-1)
                    if ((p-1)>=7 and (p-1)<=10):
                        nlist.append((p-1)+4) ###might be a three           
                    if (p-1)>10:
                        nlist.append(p-1-4)###might be a three
            if p>n:
                    nlist.append(p+1)
                    if ((p+1)>=7 and (p+1)<=10):
                        nlist.append((p+1)+4) ###might be a three           
                    if (p-1)>10:
                        nlist.append(p+1-4)###might be a three   
        Results= []
        for r in nlist:
            Results.append((r,)+self.SingleFit(image,int(r))+(imageno,))
        return max(Results,key=operator.itemgetter(1,0), default=[0,0,0,0])
    def SingleFit(self,image,n):
        self.loadDirectImageData(image)
        if n<=10: ## might be a 10
            self.closed = False
            self.makeElementsFromDeltas(np.zeros(n), 2*math.pi/9)
            self.makeTemplateFromPoints()
            optresult = self.runOptimization(x0 = np.zeros(n+1))
        else:
            self.closed = True
            self.makeElementsFromDeltas(np.zeros(n-5), 0) ###might be a three
            self.makeTemplateFromPoints()
            optresult = self.runOptimization(x0 = np.zeros(n-5+1))
            optresult.x = np.append(optresult.x,-sum(optresult.x[1:]))###might be a three
        CostFunct = optresult.fun 
        AngleResult= optresult.x 
        return (1/CostFunct,AngleResult)

def ImageFit(imageNo):
    global SaveFolderPath
    myRingFit = SAS6(nativeAngle = 2*math.pi/9,nativeLength = 4.2)
    myRingFit.loadSegmentation()  
    myRingFit.SetSavePath(SaveFolderPath)
    return myRingFit.processSingleImage(imageNo)




class MyBounds(object):
    def __init__(self, angleBounds=[-100, 100], deltaBounds = [-math.pi/10,math.pi/10]):
        self.angleBounds = angleBounds
        self.deltaBounds = deltaBounds
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        angle = x[0]
        deltas = x[1:]
    #   print(angle)
    #   print(deltas)
        tmax = bool(np.all(deltas <= self.deltaBounds[1]))
        tmin = bool(np.all(deltas >= self.deltaBounds[0]))
        amax = bool(angle <= self.angleBounds[1])
        amin = bool(angle >= self.angleBounds[0])
    #    print('{0} {1} {2} {3}'.format(tmax, tmin, amax, amin))
        return tmax and tmin and amax and amin

#myRingFit = SAS6(nativeAngle = 2*math.pi/9)

#myRingFit.makeElementsFromDeltas(np.zeros(4), 0.4)
#myRingFit.makeTemplateFromPoints()
#print(myRingFit.elementPos)
#print(np.min(myRingFit.elementPos,axis=0))
#print(np.max(myRingFit.elementPos,axis=0)-np.min(myRingFit.elementPos,axis=0))
#myRingFit.runOptimization(x0 = np.zeros(5))
#print(myRingFit.matchTemplate())
#myRingFit.showImageData()
#myRingFit.overlay_image_mask()
#myRingFit.showTemplateOverlay()
class MyTakeStep(object):
   def __init__(self, stepsize=0.5):
       self.stepsize = stepsize
   def __call__(self, x):
       s = self.stepsize
       x[0] += np.random.uniform(-10.*s, 10.*s)
       x[1:] += np.random.uniform(-s, s, x[1:].shape)
       return x
  

class ImageContainer():
    class __ImageContainer:
        def __init__(self, fileName):
            self.imageFile = h5py.File(fileName, 'r') 
            self.calculateClassAndTime()
            
        def __str__(self):
            return repr(self) + self.val
        
        def calculateClassAndTime(self):
            classes=list()
            for i in range(self.imageFile['table'].shape[0]): 
                classes.append(self.imageFile['table'][i][4].decode())     
            uniqueclasses=np.ndarray.tolist(np.unique(classes))
            #uniqueclasses = [i.decode() for i in uniqueclasses] 
            uniqueclasses = natsorted(uniqueclasses)
            #classnumber=len(uniqueclasses)
            classnum=np.zeros(len(classes))
            times= np.zeros(self.imageFile['table'].shape[0])
            for i in range(self.imageFile['table'].shape[0]): 
                times[i]=self.imageFile['table'][i][1] 
            for i, word in enumerate(classes):
                classnum[i] = uniqueclasses.index(word)+1
            self.classnum = classnum
            self.times = times
            return True
    
        def getClasses(self):
            return self.classnum
        def getTimes(self):
            return self.times
        def getImage(self,imageno):
            image1=np.asarray(self.imageFile['images'][natsorted(list(self.imageFile['images'].keys()))[imageno]]['raw'])
            image2=np.asarray(self.imageFile['images'][natsorted(list(self.imageFile['images'].keys()))[imageno]]['labeling'])#'raw'/'labeling'
            imageData  = image1*(image2/image2.max())
            #imageData = imageData+ np.min(imageData[np.nonzero(imageData)])*(1-image2/image2.max()) 
            imageData  = imageData.astype(np.uint8)
            return imageData

    instance = None
    def __init__(self, arg):
        if not ImageContainer.instance:
            ImageContainer.instance = ImageContainer.__ImageContainer(arg)
        else:
            pass
            #ImageContainer.instance.val = arg
    def __getattr__(self, name):
        return getattr(self.instance, name)


filename = os.path.abspath("ExportedForFit-RegisteredRescaled.h5")  
SaveFolderPath=Path(os.path.abspath('FitResults'))
myImageContainer = ImageContainer(filename)
results = []
nstart = 1
ntot = 11000
imageList = range(nstart,nstart+ntot)        
num_cores = 20 #multiprocessing.cpu_count()    

if __name__ == '__main__':    
    
    myPool = Pool(20) #20
    results = myPool.map(ImageFit, imageList)
    



    
#results = [ImageFit(imageNo) for imageNo in imageList ]
c=['Class','Time','PredClass','CostFunc','Angles','PartNo']
classes=myImageContainer.getClasses()
times = myImageContainer.getTimes()
df = pd.DataFrame(list(zip(classes[nstart:nstart+ntot],times[nstart:nstart+ntot],[num[0] for num in results],[num[1] for num in results],[num[2] for num in results],[num[3] for num in results])), columns=c)
df.to_csv("out.csv")
#myRingFit = SAS6(nativeAngle = math.pi/9)
#myRingFit.loadImageData('Testcase1.png')
#myRingFit.makeElementsFromDeltas(np.zeros(9))

#myRingFit.closed = True

#myRingFit.makeTemplateFromPoints()
#myRingFit.overlay_image_mask()
