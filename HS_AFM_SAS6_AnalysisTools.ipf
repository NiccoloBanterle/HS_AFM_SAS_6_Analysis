#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// FUNCTIONS FOR OPENING CLOSING (Figure 2)
Function FitOnSequence(croppedin,tres,minrad,radstep,radnum)
wave croppedin
variable tres,,minrad,radstep,radnum
make/O/N=(dimsize(croppedin,2),20) Profiles, ProfX,ProfY
make/O/N=(dimsize(croppedin,2)) Radi
variable i
for(i=0;i<dimsize(croppedin,2);i+=1)
	make/O/N=(dimsize(croppedin,0),dimsize(croppedin,1)) SigFitPlane=croppedin[p][q][i]
	RingFit(SigFitPlane,minrad,radstep,radnum,tres)
	nvar Maxradius
	wave W_ImageLineProfile,W_LineProfileX,W_LineProfileY
	Profiles[i][]=W_ImageLineProfile[q]
	ProfX[i][]=W_LineProfileX[q]
	ProfY[i][]=W_LineProfileY[q]
	Radi[i]=Maxradius
	
	print i
endfor
End

function RingFit(imagein,minrad,radstep,radnum,tres)
wave imagein
variable minrad,radstep,radnum
variable tres
variable i,rad,sump=0,currsum,maxr
maxr=Nan
for(i=0;i<radnum;i+=1)
	rad=minrad+i*radstep
	FindDxAndDyAndShift(imagein,rad,tres)
	wave xpos_sh,ypos_sh
	currsum=sumprofile(imagein,xpos_sh,ypos_sh,tres)
	if(sump<currsum)
		sump=currsum
		maxr=rad	
	endif
endfor
FindDxAndDyAndShift(imagein,maxr,tres)
wave xpos_sh,ypos_sh
variable/G Maxradius=maxr
imagelineprofile/V xwave=xpos_sh,ywave=ypos_sh, srcwave=imagein,width=1
return maxr
end

Function FindDxAndDyAndShift(filtered,radius,tres)
wave filtered 
variable radius,tres
MakeRadiusPos(radius,20)
wave xpos,ypos
make/O/N=(dimsize(xpos,0)) xpos_sh
make/O/N=(dimsize(ypos,0)) ypos_sh
make/o/n=(dimsize(filtered,0),dimsize(filtered,1)) ShiftTable=0
variable i,j
make/O/N=(dimsize(filtered,0),dimsize(xpos,0)) xpos_s
make/O/N=(dimsize(filtered,1),dimsize(ypos,0)) ypos_s
xpos_s[][]=xpos[q]+p
ypos_s[][]=ypos[q]+p
multithread ShiftTable[][]=SumProfile(filtered,row(xpos_s,p),row(ypos_s,q),tres)
imagestats ShiftTable
xpos_sh=xpos+V_maxRowLoc
ypos_sh=ypos+V_maxColLoc
variable/G Xshift=V_maxRowLoc
variable/G Yshift=V_maxColLoc
end

threadsafe function SumProfile(image,wavex,wavey,tres)
wave image, wavex,wavey
variable tres
variable inc
variable sumprof
DFREF dfSav= GetDataFolderDFR()
imagelineprofile/SC xwave=wavex,ywave=wavey, srcwave=image,width=1
wave W_ImageLineProfile
duplicate/o/free W_ImageLineProfile,W_ImageLineProfile_w
W_ImageLineProfile_w = W_ImageLineProfile>tres ? W_ImageLineProfile[p]  : Nan
wavetransform/O zapNaNs, W_ImageLineProfile_w
//imagelineprofile/V xwave=wavex,ywave=wavey, srcwave=image,width=1
//wave W_ImageLineProfile
//W_ImageLineProfile_w+=W_ImageLineProfile
if(dimsize(W_ImageLineProfile_w,0)>0)
//variable pointovertrasold
	sumprof=sum(W_ImageLineProfile_w)-(dimsize(image,0))*abs((mean(W_ImageLineProfile_w)-wavemax(image)))//-variance(W_ImageLineProfile_w)
else
	sumprof=0
endif	
SetDataFolder dfSav
return sumprof
end

function MakeRadiusPos(radius,npoint)
variable radius, npoint
make/N=(npoint)/o xpos,ypos
xpos[]=radius*cos(p*2*pi/(dimsize(xpos,0)-1))
ypos[]=radius*sin(p*2*pi/(dimsize(xpos,0)-1))
xpos[npoint-1]=radius*cos(0)
ypos[npoint-1]=radius*sin(0)
end


threadsafe function/wave row(w,n)
  wave w;
  variable n;

  make/n=0/free wc;
  if(n < dimsize(w,0))
    matrixop/free wc = row(w,n);
  endif
  redimension/N=(dimsize(w,1)) wc
  return wc;
end

Function BatchProcessRings(tres,minrad,radstep,radnum)
variable tres,minrad,radstep,radnum
variable i
string ringlist,ringinname
ringlist=Wavelist("Ring*",";","")
for(i=0;i<itemsinlist(ringlist);i+=1)
	ringinname=stringFromList(i,ringlist)
	wave ringin=$ringinname
	 FitOnSequence(ringin,tres,minrad,radstep,radnum)
	 wave Profiles, ProfX, ProfY, Radi
	 string newprofname,newxname,newyname,newradname
	 newprofname=ringinname+"Prof"
	 newxname=ringinname+"XPos"
	 newyname=ringinname+"YPos"
	 newradname=ringinname+"Rad"
	 rename Profiles, $newprofname
	 rename ProfX, $newxname
	 rename ProfY, $newyname
	 rename Radi, $newradname
	 RadiAndProfilesToTransitions($newradname,$newprofname,tres,ringinname)
	 MakeExpandedPositions($newxname,$newyname,ringinname,10)
	 AvClass(ringinname)
	 string openclosedname="openclosed"+ringinname
	 FindTransitions($newradname,$openclosedname,$ringinname)
endfor
end

Function BatchProcessExperiment()
variable noexperiments
string datafolderlist, SelectedList
string homedata=GetDataFolder(1)
datafolderlist =datafolderdir(1)
datafolderlist=datafolderlist[8,strlen(datafolderlist)-3]
SelectedList=greplist(datafolderlist,"Exp*",0,",")
NoExperiments= itemsinlist(selectedlist,",")
variable i
for(i=0;i<NoExperiments;i+=1)
	string currentexp
	currentexp=stringfromlist(i,SelectedList,",")
	setdatafolder $homedata
	setdatafolder $currentexp
	BatchProcessRings(150,5,1,5)
endfor	
end

Function RadiAndProfilesToTransitions(radius,profile,tres,basename)
wave radius,profile
variable tres
string basename
variable i,j
make/O/FREE/N=(dimsize(profile,1)) TempProf
string openclosedname="OpenClosed"+basename
make/O/n=(dimsize(radius,0)) $openclosedname
wave openclosed=$openclosedname
	for(j=0;j<dimsize(radius,0);j+=1)//for each time point
		tempprof[]=profile[j][p]
		if(OpenOrClose(tempprof,tres)==0)
			OpenClosed[j]=0
		else
			OpenClosed[j]=1
		endif	
	endfor
end

function OpenOrClose(wavein,tres)
wave wavein
variable tres
variable i,temp
temp=0
for(i=0;i<dimsize(wavein,0);i+=1)
	if(wavein[i]<=tres)
		temp=1
	endif
endfor
return temp
end



Function FindTransitions(radius,openclosed,ring)
wave radius,openclosed,ring
string basename=nameofwave(ring)
variable j
variable status,isfirst,transstarttime,transstartstatus,transstop,COCCount,OCOCount,k,maxlenghtCOC=1,maxlenghtOCO=1
	COCCount=0
	OCOCount=0
	status=openclosed[0]
	isfirst=1
	for(j=0;j<dimsize(radius,0);j+=1)
		if(status==openclosed[j])//if there is no transition do nothing
			
		else //if there is a transition
			if(isfirst==1) //check if it's the first
					isfirst=0 //change and it is not the firstanymore
					status=openclosed[j]
					transstarttime=j//start the transition
				else
					transstop=j//end of the transition started at transstarttime
					if(status==1)//it is a COC transition
						if(radius[transstarttime-1]==radius[transstop] && numtype(radius[transstop])==0 )
							COCCount+=1//count the valid COC transition 
							if((transstop-transstarttime)>maxlenghtCOC)
								maxlenghtCOC=(transstop-transstarttime)
							endif
						endif
					else //it is a OCO
						transstartstatus=0
						for(k=(transstarttime);k<(j-1);k+=1)//Check if radius between C and C stayed the same
							if(radius[k]!=radius[k+1] && numtype(radius[k])==0)
								transstartstatus=1
							endif
						endfor	
						if(transstartstatus==0)
								OCOCount+=1
								if((transstop-transstarttime)>maxlenghtOCO)
									maxlenghtOCO=(transstop-transstarttime)
								endif
						endif
					endif
					transstarttime=j
					status=openclosed[j]//change the status
			endif
		endif
	endfor
	
	//Count as Radius_Closing 
	//Save CroppedRegion COC 
//find Open->Closed->Open 
	//Check if radius Closed stayed constant 
	//Count as Radius_Opening
	//SavedCroppedRegion OCO
string COCStarttimename,COCStoptimename,COCRadiusname
string OCOStarttimename,OCOStoptimename,OCORadiusname
COCStarttimename="COCStarttime"+basename 
COCStoptimename="COCStoptime"+basename 
COCRadiusname="COCRadius"+basename 
OCOStarttimename="OCOStarttime"+basename 
OCOStoptimename="OCOStoptime"+basename 
OCORadiusname="OCORadius"+basename 
string ocoNoSpokesname,cocNosPokesname
ocoNoSpokesname="OcoNoSpokes"+basename
cocNosPokesname="CoCNoSpokes"+basename
make/O/N=(OCOCount) $ocoNoSpokesname
make/O/N=(COCCount) $cocNosPokesname
wave ocomedangle=$ocoNoSpokesname
wave cocmedangle=$cocNosPokesname

make/O/N=(COCCount) $COCStarttimename=Nan,$COCStoptimename=Nan,$COCRadiusname=Nan
make/O/N=(OCOCount) $OCOStarttimename=Nan,$OCOStoptimename=Nan,$OCORadiusname=Nan

wave COCStarttime=$COCStarttimename
wave COCStoptime=$COCStoptimename
wave COCRadius=$COCRadiusname
wave OCOStarttime=$OCOStarttimename
wave OCOStoptime=$OCOStoptimename
wave OCORadius=$OCORadiusname
string classname=basename+"class"
wave class=$classname
COCCount=0
OCOCount=0
	status=openclosed[0]
	isfirst=1
	for(j=0;j<dimsize(radius,0);j+=1)
		if(status==openclosed[j])//if there is no transition do nothing
			
		else //if there is a transition
			if(isfirst==1) //check if it's the first
					isfirst=0
					status=openclosed[j]//change and it is not the firstanymore
					transstarttime=j
				else
					transstop=j//end of the transition started at transstarttime
					if(status==1)//it is a COC transition
						if(radius[transstarttime-1]==radius[transstop] && numtype(radius[transstop])==0)
							COCStarttime[COCCount]=transstarttime
							COCStoptime[COCCount]=transstop
							COCRadius[COCCount]=radius[transstop]
							cocmedangle[COCCount]=class[transstarttime]
							COCCount+=1 //count the valid COC transition 
						endif
					else //it is a OCO
						transstartstatus=0
						for(k=(transstarttime);k<(j-1);k+=1)//Check if radius between C and C stayed the same
							if(radius[k]!=radius[k+1] && numtype(radius[j])==0)
								transstartstatus=1
							endif
						endfor	
						if(transstartstatus==0)
								OCOStarttime[OCOCount]=transstarttime
								OCOStoptime[OCOCount]=transstop
								OCORadius[OCOCount]=radius[transstop-1]
								ocomedangle[OCOCount]=class[transstarttime+1]
								OCOCount+=1
						endif
					endif
					transstarttime=j
					status=openclosed[j]//change the status
			endif
		endif
	endfor
string ocolenghtname,coclenghtname
ocolenghtname="OCOLenght"+basename
coclenghtname="COCLenght"+basename

make/O/N=(dimsize(OCOStarttime,0)) $OCOLenghtname=OCOStoptime-OCOStarttime
make/O/N=(dimsize(COCStarttime,0)) $COCLenghtname=COCStoptime-COCStarttime



end

function RefineTransitionSimmetry(ring,addition)
	wave ring
	variable addition
	string ringname
	string nameofcocstart,nameofcocstop, nameofocostart,nameofocostop
	wave cocStart=$nameofcocstart
	wave cocStop=$nameofcocstop
	wave ocostart=$nameofocostart
	wave ocostop=$nameofocostop
	string profilexname,profileyname
	string nameofprofilex,nameofprofiley
	nameofprofilex=nameofwave(ring)+"XPos"
	nameofprofiley=nameofwave(ring)+"YPos"
	wave Profilex=$nameofprofilex
	wave ProfileY=$nameofprofiley
	MakeExpandedPositions(Profilex,ProfileY,nameofwave(ring),addition)
end

function MakeExpandedPositions(XPosall,YPosall,basename,addition)
wave xposall,yposall
string basename
variable addition 
	variable i
	make/FREE/N=(dimsize(xposall,1)) xpost,ypost
	string nameofx,nameofy,nameofprof	, nameofclassouter,mediananglename,allpeakinfoname
	nameofx=nameofwave(XPosall)+"Outer"
	nameofy=nameofwave(YPosall)+"Outer"
	nameofclassouter=basename+"Spokesno"
	nameofprof=basename+"ProfOuter"
	mediananglename=basename+"MediAngle"
	allpeakinfoname=basename+"AllPeaks"
	make/O/N=(dimsize(xposall,0),4*dimsize(xposall,1)) $nameofx,$nameofy,$nameofprof
	make/O/N=(dimsize(xposall,0)) $nameofclassouter, $mediananglename
	make/O/N=(11,5,dimsize(xposall,0)) $allpeakinfoname
	wave newx=$nameofx
	wave newy=$nameofy
	wave newprof=$nameofprof
	wave rings=$basename
	wave spokesno=$nameofclassouter
	wave medianangle=$mediananglename
	wave allpeaks=$allpeakinfoname
	make/O/N=(dimsize(rings,0),dimsize(rings,1))/free tempring
	for(i=0;i<dimSize(xposall,0);i+=1)
		xpost=xposall[i][p]
		ypost=yposall[i][p]
		variable xcenter,ycenter,radius
		xcenter=(wavemin(xpost)+wavemax(xpost))/2
		ycenter=(wavemin(ypost)+wavemax(ypost))/2
		radius=(wavemax(xpost)-xcenter)
		MakeRadiusPos(radius+addition,4*dimsize(xposall,1))
		wave xpos,ypos
		newx[i][]=xpos[q]+xcenter
		newy[i][]=ypos[q]+ycenter
		xpos+=xcenter
		ypos+=ycenter
		tempring[][]=rings[p][q][i]
		imagelineprofile xwave=xpos,ywave=ypos, srcwave=tempring,width=2
		wave W_ImageLineProfile
		newprof[i][]=W_ImageLineProfile[q]
		spokesno[i]=AutoFindPeaks(W_ImageLineProfile,0,dimsize(W_ImageLineProfile,0)-1,0.1,2,11)
		wave W_AutoPeakInfo
		allpeaks[0,dimsize(W_AutoPeakInfo,0)-1][][i]=W_AutoPeakInfo[p][q]
		make/FREE/N=(dimsize(W_AutoPeakinfo,0)) PeakPos,PeakDist
		PeakPos[]=W_AutoPeakInfo[p][0]
		sort PeakPos,Peakpos
		PeakDist[0,dimsize(peakpos,0)-2]=abs(PeakPos[p]-peakpos[p+1])
		PeakDist[dimsize(peakpos,0)-1]=abs(peakpos[0])+abs(dimsize(W_ImageLineProfile,0)-1-peakpos[dimsize(peakpos,0)-1])
		medianangle[i]=(median(PeakDist)/dimsize(W_ImageLineProfile,0))*360
endfor
end

function AvClass(basename)
string basename
string openclosedname,mediananglename,classname
openclosedname="openclosed"+basename
mediananglename=basename+"mediangle"
classname=basename+"class"
wave opclosed=$openclosedname
wave medangle=$mediananglename
make/O/N=(dimsize(opclosed,0)) $classname
wave class=$classname
variable count=0
variable i
make/FREE/N=(dimsize(opclosed,0)) temptrans
for(i=1;i<dimsize(opclosed,0);i+=1)
	if(opclosed[i]!=opclosed[i-1])
		count+=1
	endif	
	temptrans[i]=count
endfor

for(i=0;i<count;i+=1)
make/FREE/N=(dimsize(opclosed,0)) tempmedian
	tempmedian = temptrans[p]==i ? medangle[p] : Nan
	wavetransform/o zapnans, tempmedian
	class = temptrans[p]==i ? median(tempmedian) : class[p]
	class = round(360/class[p])
endfor
end

// FUNCTIONS FOR KINETICS (Figure 3)

#include <KBColorizeTraces>
#include <Global Fit 2>
function ConvertProbToClass(timewave,probmatrix)
wave timewave, probmatrix
make/O/N=(wavemax(timewave)+1,dimsize(probmatrix,1)) TimeMatrix=0
variable,class,i
for(class=0;class<dimsize(probmatrix,1);class+=1)
	for(i=0;i<(dimsize(probmatrix,0)-1);i+=1)
		TimeMatrix[timewave[i]][class]+=probmatrix[i][class]
	endfor
endfor
end

//Analysis
//ProbTable[][0]=table_ProbabilityofLabel1[p]
//•ProbTable[][1]=table_ProbabilityofLabel2[p]
//•ProbTable[][2]=table_ProbabilityofLabel3[p]
//•ProbTable[][3]=table_ProbabilityofLabel4[p]
//•ProbTable[][4]=table_ProbabilityofLabel5[p]
//•ProbTable[][5]=table_ProbabilityofLabel6[p]
//•ProbTable[][6]=table_ProbabilityofLabel7[p]
//•ProbTable[][7]=table_ProbabilityofLabel8[p]
//•ProbTable[][8]=table_ProbabilityofLabel9[p]
//•ProbTable[][9]=table_ProbabilityofLabel10[p]
//•ProbTable[][10]=table_Probabilityof7closed[p]
//•ProbTable[][11]=table_Probabilityof8closed[p]
//•ProbTable[][12]=table_Probabilityof9closed[p]
//•ProbTable[][13]=table_Probabilityof10closed[p]
//•ConvertProbToClass(table_timestep,ProbTable)

function GenerateError(timewave,probmatrix)
wave timewave, probmatrix
make/Free/O/N=(wavemax(timewave)+1,dimsize(probmatrix,1)) TempErrorMatrix=0
make/O/N=(wavemax(timewave)+1,dimsize(probmatrix,1)) ErrorMatrix=0
make/O/N=(wavemax(timewave)+1) Objectnumber=0
variable,class,i
for(class=0;class<dimsize(probmatrix,1);class+=1)
	for(i=0;i<(dimsize(probmatrix,0)-1);i+=1)
		TempErrorMatrix[timewave[i]][class]+=probmatrix[i][class]
		Objectnumber[timewave[i]]+=1
	endfor
endfor
Objectnumber/=dimsize(probmatrix,1)
ErrorMatrix[][]=sqrt(TempErrorMatrix[p][q]*(Objectnumber[p]-TempErrorMatrix[p][q])/Objectnumber[p])
end

function GenerateIndividualWaves(matrixin)
wave matrixin
variable i
string nameofnew
for(i=0;i<dimsize(matrixin,1);i+=1)
	nameofnew="Conc"+num2str(i+1)
	make/O/N=(dimsize(matrixin,0)) $nameofnew
	wave current=$nameofnew
	current[]=matrixin[p][i]
endfor
end

function GenerateIndividualErrors(matrixin)
wave matrixin
variable i
string nameofnew
for(i=0;i<dimsize(matrixin,1);i+=1)
	nameofnew="Error"+num2str(i+1)
	make/O/N=(dimsize(matrixin,0)) $nameofnew
	wave current=$nameofnew
	current[]=matrixin[p][i]
endfor
end

Function TotalConcentration(timematrixin)
wave timematrixin
make/O/N=(dimsize(timematrixin,0)) TotConc=0
TotConc[]=timematrixin[p][0]+2*timematrixin[p][1]+3*timematrixin[p][2]+4*timematrixin[p][3]+5*timematrixin[p][4]+6*timematrixin[p][5]+7*timematrixin[p][6]+8*timematrixin[p][7]+9*timematrixin[p][8]+10*timematrixin[p][9]
TotConc[]+=7*timematrixin[p][10]+8*timematrixin[p][11]+9*timematrixin[p][12]+10*timematrixin[p][13]
end

Function FitRingAssemblyComplete78910(pw, yw,xw) : FitFunc
	Wave pw		// pw[0] = kon, pw[1] = koff, pw[2] = kringon8, pw[3]=kingoff8, p[4]=SolConc,p[5]=t0,p[6]=DiffSteep,p[7]=AdsMax,p[8]=Des, p[9]=AdRate, p[10]=Kringon9,p[11]=kringoff9, p[12]=Kringon10, p[13]=kringoff10, p[14]=kringoff7,,p[15]=kringoff7, p[16]=TimePoints, p[17]=StateNo, p[18]=OligomerState
	wave yw		
	wave xw
	make/O/N=(16) KineticConstantsVC=PW[p] 
	make/O/N=(PW[16],PW[17]) FullmatrixVC=0
 	IntegrateODE /M=0 RingAssemblyKinVarConc78910, KineticConstantsVC, FullMatrixVC
 	yw[]=FullMatrixVC[p][pw[18]]
 	return 0
End

Function RingAssemblyKinVarConc78910(pw, tt, yw, dydt)
	Wave pw						// pw[0] = kon, pw[1] = koff, pw[2] = kringon8, pw[3]=kingoff8, p[4]=SolConc,p[5]=t0,p[6]=DiffSteep,p[7]=AdsMax,p[8]=Des, p[9]=AdRate, p[10]=Kringon9,p[11]=kringoff9, p[12]=Kringon10, p[13]=kringoff10 p[14]=Kringon7, p[15]=kringoff7
	Variable tt					// time value at which to calculate derivatives
	Wave yw						// yw[0]-yw[8] containing concentrations of C1-C9
	Wave dydt						// wave to receive dA/dt, dB/dt etc. (output)
	variable j,l,sum1,sum2,sum3,sum4,sum5,Csol,Csurf
	
Csol=pw[4]*(Erf(pw[6]*(tt-pw[5]))+1)
Csurf=0
variable i
for(i=0;i<dimsize(yw,0);i+=1)
	CSurf+=yw[i]*(i<=10 ? i+1:i+1-4) //multiplies for the 
endfor
	for(j=0;j<=9;j+=1)
		sum1=0
		for(l=0;l<=(j-1);l+=1)
			sum1+=yw[l]*yw[j-l-1]
		endfor
		sum2=0
		for(l=j+1;l<=9;l+=1)
			sum2+=yw[l]
		endfor
		sum3=0		
		for(l=0;l<=(9-j-1);l+=1)
			sum3+=yw[l]
		endfor	
		dydt[j]= pw[0]*sum1+2*pw[1]*sum2-pw[1]*(j)*yw[j]-pw[0]*2*yw[j]*sum3-(j==7)*pw[2]*yw[7]+pw[3]*yw[11]*(j==7)-(j==8)*pw[10]*yw[8]+pw[11]*yw[12]*(j==8)-(j==9)*pw[12]*yw[9]+(j==9)*pw[13]*yw[13]+(j==0)*CSol*pw[9]*(pw[7]-Csurf)-pw[8]*Csurf*(j==0)-(j==6)*pw[14]*yw[6]+(j==6)*pw[15]*yw[10]
	endfor
	dydt[10]=pw[14]*yw[6]-pw[15]*yw[10]//7closed
	dydt[11]=pw[2]*yw[7]-pw[3]*yw[11]//8closed
	dydt[12]=pw[10]*yw[8]-pw[11]*yw[12]//9closed
	dydt[13]=pw[12]*yw[9]-pw[13]*yw[13]//9closed
	return 0
End

function FitOverallConcComplete78910(pw, yw,xw): FitFunc
Wave pw		// pw[0] = kon, pw[1] = koff, pw[2] = kringon8, pw[3]=kingoff8, p[4]=SolConc,p[5]=t0,p[6]=DiffSteep,p[7]=AdsMax,p[8]=Des,
// p[9]=AdRate, p[10]=Kringon9,p[11]=kringoff9, p[12]=Kringon10, p[13]=kringoff10, p[14]=kringoff7,,p[15]=kringoff7, p[16]=TimePoints, p[17]=StateNo, p[18]=OligomerState
	wave yw		
	wave xw
 	make/O/N=(19) KineticConstantsTotConcVC=PW[p] 
	make/O/N=(PW[16],PW[17]) FullmatrixTC=0, WeightFMVC=0 
 	IntegrateODE /M=0 RingAssemblyKinVarConc78910, KineticConstantsTotConcVC, FullmatrixTC
 	WeightFMVC=FullmatrixTC
 	WeightFMVC[][0,9]*=q+1
 	WeightFMVC[][10]*=7
 	WeightFMVC[][11]*=8
 	WeightFMVC[][12]*=9
 	WeightFMVC[][13]*=10
 	sumdimension/DEST=yw/D=1 WeightFMVC 	
 	return 0
End

#include <ImageSlider>
#include <Multi-peak fitting 2.0>
#include <Peak AutoFind>
#include <Scatter Plot Matrix 2>
#include <TransformAxis1.2>



// FUNCTIONS FOR PHASE DIAGRAMS (Figure 4)


function SimulateRingAssemblyComplete(inputparam)//needed
wave inputparam
make/O/N=(16) SimKineticConstantsVC=inputparam[p] 
make/O/N=(inputparam[16],inputparam[17]) FullmatrixSim=0
IntegrateODE/M=0 RingAssemblyKinVarConc78910, SimKineticConstantsVC, FullmatrixSim
end 

function SimulateRingAssemblyCompleteConsConc(inputparam)//needed
wave inputparam
make/O/N=(11) SimKineticConstantsVC=inputparam[p] 
make/O/N=(inputparam[11],inputparam[12]) FullmatrixSim=0
FullMatrixSim[0][0]=inputparam[4]
IntegrateODE/M=0 RingAssemblyKinConstConc78910, SimKineticConstantsVC, FullmatrixSim
end 

Function PhaseDiagram(paramringassembly,param1,param2,p1st,p1bin,p1num,p2st,p2bin,p2num)
wave paramringassembly
variable param1,param2,p1st,p1bin,p1num,p2st,p2bin,p2num
duplicate/O paramringassembly,  paramringassemblyvar
variable i,j
make/O/N=(p1num,p2num) w_PhaseDiagramAv,w_PhaseDiagramStd,w_RingEfficiency,w_9RingEfficiency,w_FinalConc
for(i=0;i<p1num;i+=1)
	for(j=0;j<p2num;j+=1)
		paramringassemblyvar[param1]=p1st+p1bin*i
		paramringassemblyvar[param2]=p2st+p2bin*j
		SimulateRingAssemblyComplete(paramringassemblyvar)
		wave fullmatrixsim 
		w_PhaseDiagramAv[i][j]=(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]*8+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]*9+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]*10)/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12])
		w_PhaseDiagramStd[i][j]=((fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]*(8-w_PhaseDiagramAv[i][j])^2+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]*(9-w_PhaseDiagramAv[i][j])^2+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]*(10-w_PhaseDiagramAv[i][j])^2))/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12])
		w_RingEfficiency[i][j]=(7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][13])/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+2*fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+3*fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+4*fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+5*fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+6*fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][13])
		w_9RingEfficiency[i][j]=(fullmatrixsim[dimsize(fullmatrixsim,0)-1][11])/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12])
		w_FinalConc[i][j]=fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+2*fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+3*fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+4*fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+5*fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+6*fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]
	endfor
endfor
end



Function RingAssemblyKinConstConc78910(pw, tt, yw, dydt)//needed
	Wave pw						// pw[0] = kon, pw[1] = koff, pw[2] = kringon8, pw[3]=kingoff8, p[4]=Conc, p[5]=Kringon9,p[6]=kringoff9, p[7]=Kringon10, p[8]=kringoff10 p[9]=Kringon7, p[10]=kringoff7
	Variable tt					// time value at which to calculate derivatives
	Wave yw						// yw[0]-yw[9] containing concentrations of C1-C10  
	Wave dydt						// wave to receive dA/dt, dB/dt etc. (output)
	variable j,l,sum1,sum2,sum3,sum4,sum5,Csurf
	variable i
	for(j=0;j<=9;j+=1)
		sum1=0
		for(l=0;l<=(j-1);l+=1)
			sum1+=yw[l]*yw[j-l-1]
		endfor
		sum2=0
		for(l=j+1;l<=9;l+=1)
			sum2+=yw[l]
		endfor
		sum3=0		
		for(l=0;l<=(9-j-1);l+=1)
			sum3+=yw[l]
		endfor	
		dydt[j]= pw[0]*sum1+2*pw[1]*sum2-pw[1]*(j)*yw[j]-pw[0]*2*yw[j]*sum3-(j==7)*pw[2]*yw[7]+pw[3]*yw[11]*(j==7)-(j==8)*pw[5]*yw[8]+pw[6]*yw[12]*(j==8)-(j==9)*pw[7]*yw[9]+(j==9)*pw[8]*yw[13]-(j==6)*pw[9]*yw[6]+(j==6)*pw[10]*yw[10]
	endfor
	dydt[10]=pw[9]*yw[6]-pw[10]*yw[10]//7closed
	dydt[11]=pw[2]*yw[7]-pw[3]*yw[11]//8closed
	dydt[12]=pw[5]*yw[8]-pw[6]*yw[12]//9closed
	dydt[13]=pw[7]*yw[9]-pw[8]*yw[13]//9closed
	return 0
End

Function PhaseDiagramConstConc(paramringassembly,param1,param2,p1st,p1bin,p1num,p2st,p2bin,p2num)
wave paramringassembly
variable param1,param2,p1st,p1bin,p1num,p2st,p2bin,p2num
duplicate/O paramringassembly,  paramringassemblyvar
variable i,j
make/O/N=(p1num,p2num) w_PhaseDiagramAv,w_PhaseDiagramStd,w_RingEfficiency,w_9RingEfficiency,w_FinalConc
for(i=0;i<p1num;i+=1)
	for(j=0;j<p2num;j+=1)
		paramringassemblyvar[param1]=p1st+p1bin*i
		paramringassemblyvar[param2]=p2st+p2bin*j
		SimulateRingAssemblyCompleteConsConc(paramringassemblyvar)
		wave fullmatrixsim 
		w_PhaseDiagramAv[i][j]=(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]*7+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]*8+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]*9+fullmatrixsim[dimsize(fullmatrixsim,0)-1][13]*10)/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][13])
		w_PhaseDiagramStd[i][j]=((fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]*(8-w_PhaseDiagramAv[i][j])^2+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]*(9-w_PhaseDiagramAv[i][j])^2+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]*(10-w_PhaseDiagramAv[i][j])^2))/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12])
		w_RingEfficiency[i][j]=(7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][13])/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+2*fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+3*fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+4*fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+5*fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+6*fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][13])
		w_9RingEfficiency[i][j]=(fullmatrixsim[dimsize(fullmatrixsim,0)-1][11])/(fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+fullmatrixsim[dimsize(fullmatrixsim,0)-1][12])
		w_FinalConc[i][j]=fullmatrixsim[dimsize(fullmatrixsim,0)-1][0]+2*fullmatrixsim[dimsize(fullmatrixsim,0)-1][1]+3*fullmatrixsim[dimsize(fullmatrixsim,0)-1][2]+4*fullmatrixsim[dimsize(fullmatrixsim,0)-1][3]+5*fullmatrixsim[dimsize(fullmatrixsim,0)-1][4]+6*fullmatrixsim[dimsize(fullmatrixsim,0)-1][5]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][6]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][7]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][8]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][9]+7*fullmatrixsim[dimsize(fullmatrixsim,0)-1][10]+8*fullmatrixsim[dimsize(fullmatrixsim,0)-1][11]+9*fullmatrixsim[dimsize(fullmatrixsim,0)-1][12]+10*fullmatrixsim[dimsize(fullmatrixsim,0)-1][13]
	endfor
endfor
end

Function OpeningFactorSimulation(paramringassembly,OpfactorSt,OpFactorBin,OpFactorNum,timeat)
 wave paramringassembly
 variable OpfactorSt,OpFactorBin,OpFactorNum,timeat
 	duplicate/O paramringassembly,paramringassemblyOpFac
	variable i,opfact,totclosed
	make/O/N=(OpFactorNum) No7,No8,No9,No10
	for(i=0;i<OpFactorNum;i+=1)
		opfact=OpfactorSt+i*OpFactorBin //3,6,8,10
		paramringassemblyOpFac[3]=paramringassembly[3]*opfact
		paramringassemblyOpFac[6]=paramringassembly[6]*opfact
		paramringassemblyOpFac[8]=paramringassembly[8]*opfact
		paramringassemblyOpFac[10]=paramringassembly[10]*opfact
		SimulateRingAssemblyCompleteConsConc(paramringassemblyOpFac)
		wave fullmatrixsim 
		totclosed=fullmatrixsim[timeat][10]+fullmatrixsim[timeat][11]+fullmatrixsim[timeat][12]+fullmatrixsim[timeat][13]
		No7[i]=fullmatrixsim[timeat][10]/(totclosed)
		No8[i]=fullmatrixsim[timeat][11]/(totclosed)
		No9[i]=fullmatrixsim[timeat][12]/(totclosed)
		No10[i]=fullmatrixsim[timeat][13]/(totclosed)
	endfor
end


// FUNCTIONS FOR POLYMER CHAIN FITTING (Figure 5)

function SingClassImages(alignedstack,popwave,n)
wave alignedstack,popwave
variable n
variable i
variable counter=0
for(i=0;i<dimsize(popwave,0);i+=1)
	if(popwave[i]==n)
		counter+=1
	endif
endfor
string namepopstack
namepopstack= "waveStack"+num2str(n)
make/O/N=(dimsize(alignedstack,0),dimsize(alignedstack,1),counter) $namepopstack  
wave subpop=$namepopstack 
counter=0
for(i=0;i<dimsize(popwave,0);i+=1)
	if(popwave[i]==n)
		subpop[][][counter]=alignedstack[p][q][i]
		counter+=1
	endif
endfor
imagetransform averageimage,subpop
wave M_aveimage
string newavgname
newavgname = namepopstack+ "_Avg"
make/O/N=(dimsize(alignedstack,0),dimsize(alignedstack,1)) $newavgname 
wave newavg=$newavgname
newavg= M_aveimage
imageregistration/stck/refm=0/tstm=0 testwave=subpop,refwave=newavg
wave M_regout
string imageregstackname=namepopstack+"Rereg"
killwaves/z $imageregstackname
rename M_regout, $imageregstackname
imagetransform averageimage, M_RegOut
wave M_aveimage
string newregnewavg=imageregstackname+"Avg"
killwaves/z $newregnewavg
rename  M_aveimage,$newregnewavg
end

Function RealignImages(objwave,realignedmatrix,xposmatrix,yposmatrix,Xshiftwave,YshiftWave)
wave objwave,realignedmatrix,xposmatrix,yposmatrix,Xshiftwave,YshiftWave
string nameofrealigned
nameofrealigned = nameofwave(objwave)+"_realign"
duplicate/O objwave,$nameofrealigned
wave realigned=$nameofrealigned
variable i
for(i=0;i<dimsize(objwave,2);i+=1)
	imagetransform/P=(i) getplane,objwave
	wave  M_ImagePlane
	imagetransform/IOFF={-xposmatrix[i][0]+dimsize(M_ImagePlane,0)/2,-yposmatrix[i][0]+dimsize(M_ImagePlane,1)/2,0} offsetimage, M_ImagePlane
	wave M_OffsetImage
	
	if(realignedmatrix[i][0]==0)
		duplicate/o M_OffsetImage,M_RotatedImage
	else
		imagerotate/e=0/A=((realignedmatrix[i][0]*360)/(2*pi)) M_OffsetImage
		wave M_RotatedImage
	endif
	M_ImagePlane[][]=M_RotatedImage[p+(dimsize(M_RotatedImage,0)-dimsize(M_ImagePlane,0))/2][q+(dimsize(M_RotatedImage,1)-dimsize(M_ImagePlane,1))/2]
	imagetransform/IOFF={Xshiftwave[i],yshiftwave[i],0} offsetimage, M_ImagePlane
	wave M_OffsetImage
	realigned[][][i]=0
	realigned[][][i]=M_OffsetImage[p][q]
endfor
end

function ExtractFittedImages(objectwave,imagenumwave)
wave objectwave,imagenumwave
variable i
string nameofnewwave
nameofnewwave="Sel_"+nameofwave(objectwave)
make/O/N=(dimsize(objectwave,0),dimsize(objectwave,1),dimsize(imagenumwave,0)) $nameofnewwave
wave Selwave=$nameofnewwave
variable count=0
for(i=0;i<dimsize(imagenumwave,0);i+=1)
Selwave[][][count]=objectwave[p][q][imagenumwave[i]]
count+=1
endfor
end

Function AlignRingToTemplate(anglematrix,cost,popwave,xpostemp,ypostemp,n,tres)
wave anglematrix,cost,popwave,xpostemp,ypostemp
variable n,tres
make/O/N=(dimsize(anglematrix,0),dimsize(anglematrix,1)) AngleMatrixRealigned=Anglematrix
make/O/N=(dimsize(anglematrix,0)) Xshift,Yshift

variable i,count=0
make/N=(1,n)/o CurrentPoly
make/O/N=3 ShiftAndFirstAngle
for(i=0;i<dimsize(anglematrix,0);i+=1)
	if(popwave[i]==n && cost[i]>tres)
		CurrentPoly[0][]=anglematrix[i][q]
		ShiftAndFirstAngle=0
		AlignToTemp(shiftandFirstAngle,currentpoly,xpostemp,ypostemp,n)
		wave shiftandFirstAngle
		AngleMatrixRealigned[i][0]=ShiftAndFirstAngle[0]
		Xshift[i]=ShiftAndFirstAngle[1]
		Yshift[i]=ShiftAndFirstAngle[2]
		count+=1
	endif
endfor
string NewxPosName,NewYposname
NewxPosName="RealignedX"+num2str(n)
NewyPosName="RealignedY"+num2str(n)

make/O/N=(count,n+1) $NewxPosName,$NewyPosName
wave NewXpos=$NewxPosName
wave NewYpos=$NewyPosName
count=0
for(i=0;i<dimsize(anglematrix,0);i+=1)
	if(popwave[i]==n && cost[i]>tres)
		CurrentPoly[0][]=AngleMatrixRealigned[i][q]
		PlotAngFit(currentpoly,0,3.9,2*pi/9)
		wave xpos,ypos
		xpos+=Xshift[i]
		ypos+=yshift[i]
		NewXpos[count][]=xpos[q]
		NewYpos[count][]=ypos[q]
		count+=1
	endif
endfor
end

function AlignToTemp(shiftandFirstAngle,currentpoly,xpostemp,ypostemp,n)
wave shiftandFirstAngle,currentpoly,xpostemp,ypostemp
variable n
duplicate/o xpostemp,xpostempp
duplicate/o ypostemp,ypostempp
make/O/N=1 w
w[0]=n
make/O/N=3 RX
RX[0]=100
RX[1]=100
RX[2]=100
shiftandFirstAngle=0
Optimize/D=4/I=10000/M = {3, 0 }/S=10/X=shiftandFirstAngle/Y=3/R=Rx SingleDistance,w
End

function GenerateTemplate(n,avAngle)
variable n,avangle
make/O/N=1 w=n
string nameTemp
nameTemp="Template"+num2str(n)
make/o/N=(n) $nameTemp
wave InitialAngle=$nameTemp
InitialAngle=avAngle-2*pi/9
Optimize/S =0.5/X=initialangle TotalDistance,w
make/o/N=(1,n) Dummy
Dummy[0][]=initialangle[q] 
PlotAngFit(Dummy,0,4.02,avAngle)
wave xpos,ypos
string nameTempX,nameTempY
nameTempX=nameTemp+"X"
nameTempY=nameTemp+"Y"
duplicate xpos,$nameTempX
duplicate ypos,$nameTempY
end

function GenerateTemplateRef(n)
variable n
make/O/N=1 w=n
string nameTemp
nameTemp="TemplateR"+num2str(n)
make/o/N=(n) $nameTemp
wave InitialAngle=$nameTemp
InitialAngle=0
Optimize/M={1,1}/X=initialangle TotalDistanceRef,w
make/o/N=(1,n) Dummy
Dummy[0][]=initialangle[q] 
PlotAngFit(Dummy,0,4.02,2*pi/9)
wave xpos,ypos
string nameTempX,nameTempY
nameTempX=nameTemp+"X"
nameTempY=nameTemp+"Y"
duplicate/O xpos,$nameTempX
duplicate/O ypos,$nameTempY
end

Function PlotAngFit(anglematrix,n,lsize,mainang)
wave anglematrix
variable n,lsize,mainang
make/O/N=(dimsize(anglematrix,1)) SingAngles=Nan
SingAngles[]=anglematrix[n][p]
wavetransform zapnans, singangles
make/O/N=(dimsize(singangles,0)+1) Xpos,Ypos
variable i,prevang
Xpos[0]=0
Ypos[0]=0
prevang=mainang
for(i=1;i<dimsize(singangles,0)+1;i+=1)
Xpos[i]=Xpos[i-1]+cos(prevang+singangles[i-1]+mainang)*lsize
Ypos[i]=Ypos[i-1]+sin(prevang+singangles[i-1]+mainang)*lsize
prevang+=singangles[i-1]+mainang
endfor
end

function SingleDistance(w,shiftandFirstAngle)
wave w
wave shiftandFirstAngle
variable n
n=w[0]
wave currentpoly
wave xpostempp,ypostempp
duplicate/O currentpoly,currentpolym
currentpolym[0][0]=shiftandFirstAngle[0]
PlotAngFit(currentpolym,0,3.9,2*pi/9)
wave xpos,ypos
xpos+=shiftandFirstAngle[1]
ypos+=shiftandFirstAngle[2]
variable Totdist=0,i
for(i=0;i<dimsize(xpostempp,0);i+=1)
		TotDist+=((xpos[i]-xpostempp[i])^2+(ypos[i]-ypostempp[i])^2)
endfor
return TotDist
end
// FUNCTIONS FOR 3D MODELS (Figure 6)

function GenerateChiralTurnRanVarTurnShift(nlayer,L,H,sigma,Amp,shiftst,shiftbin,shiftnum)
	variable nlayer,L,H,sigma,Amp,shiftst,shiftbin,shiftnum
	variable ncontact,n,Shift,distancetot,maxshift,maxcontact,mindistance
	make/O/N=(9,nlayer+1) ChiralTurnVarMatrixAngle=Nan
	make/free/O/N=9 Angles1,Angles2,angles2max
	variable count,shift2
	for(count=0;count<(nlayer);count+=1)
		maxcontact=0
		mindistance=10000
			for(shift=shiftst;shift<2*(shiftst+shiftbin*shiftnum);shift+=shiftbin)
				if(shift<(shiftst+shiftbin*shiftnum))
					if(count==0)
						for(n=0;n<9;n+=1)
								Angles1[n]= gnoise(sigma)+Amp*cos(1.45*(mod((n+1),9))-0.73)	
								Angles2[n]= gnoise(sigma)+Amp*cos(1.45*(mod((n+1)+shift,9))-0.73)
								ChiralTurnvarMatrixAngle[][0]=Angles1[p]
								//ChiralTurnvarMatrixAngle[][1]=Angles2[p]
						endfor
					else
						angles1=ChiralTurnvarMatrixAngle[p][count]
						for(n=0;n<9;n+=1)
							Angles2[n]= gnoise(sigma)+Amp*cos((1.45*(mod((n+1)+shift,9))-0.73))
						endfor
					endif
					ncontact=0
					distancetot=0
					for(n=0;n<9;n+=1)
						ncontact+=IsOverlapping(H,L,angles1[n],angles2[n])
						distancetot+=Distance(H,L,angles1[n],angles2[n])
					endfor
					if(ncontact>=maxcontact  )//maxcontact)//&& (distancetot)<mindistance
						if(ncontact>maxcontact)
							maxcontact=ncontact
							maxshift=shift
							mindistance=distancetot
							ChiralTurnvarMatrixAngle[][count+1]=Angles2[p]
						else
							if((distancetot)<mindistance)
								maxcontact=ncontact
								mindistance=distancetot
								maxshift=shift
								ChiralTurnvarMatrixAngle[][count+1]=Angles2[p]
							endif
						endif
					endif
				else
					shift2=shift-(shiftst+shiftbin*shiftnum)
					if(count==0)
						for(n=0;n<9;n+=1)
								Angles1[n]= gnoise(sigma)+Amp*cos(1.45*(mod((n+1),9))-0.73)	
								Angles2[n]= gnoise(sigma)-Amp*cos(1.45*(mod((n+1)+shift2,9))-0.73)
								ChiralTurnvarMatrixAngle[][0]=Angles1[p]
								//ChiralTurnvarMatrixAngle[][1]=Angles2[p]
						endfor
					else
						angles1=ChiralTurnvarMatrixAngle[p][count]
						for(n=0;n<9;n+=1)
							Angles2[n]= gnoise(sigma)-Amp*cos((1.45*(mod((n+1)+shift2,9))-0.73))
						endfor
					endif
					ncontact=0
					distancetot=0
					for(n=0;n<9;n+=1)
						ncontact+=IsOverlapping(H,L,angles1[n],angles2[n])
						distancetot+=Distance(H,L,angles1[n],angles2[n])
					endfor
					if(ncontact>=maxcontact  )//maxcontact)//&& (distancetot)<mindistance
						if(ncontact>maxcontact)
							maxcontact=ncontact
							maxshift=shift
							mindistance=distancetot
							ChiralTurnvarMatrixAngle[][count+1]=Angles2[p]
						else
							if((distancetot)<mindistance)
								maxcontact=ncontact
								mindistance=distancetot
								maxshift=shift
								ChiralTurnvarMatrixAngle[][count+1]=Angles2[p]
							endif
						endif
					endif
			endif	
			endfor
		endfor	
end

Function IsOverlapping(H,L,a1,b1)
variable H,L,a1,b1
if(sign(a1)==sign(b1) || a1==b1)
	return 0	
else
	variable Xint
	Xint=H/(tan(a1)-tan(b1))
	if(Xint<0 || Xint>min(L*cos(a1),L*cos(b1)))
		return 0	
	else
		return 1
	endif		
endif
end

Function Distance(H,L,a1,b1)
variable H,L,a1,b1

	variable x1,x2,y1,y2
	x1=L*cos(a1)
	x2=L*cos(b1)
	y1=L*sin(a1)
	y2=H+L*sin(b1)
		return sqrt((x1-x2)^2+(y1-y2)^2)

end

Function UnimodalRandom(L,H,minspokes,sigma)
	variable L,H,minspokes,sigma
	variable count=0,ncontact,n
	make/O/N=(9,400) UnimodalMatrixAngle=Nan
	do
		make/O/N=9 Angles1,Angles2
		if(count==0)
			Angles1=gnoise(sigma)
			Angles2=gnoise(sigma)
			UnimodalMatrixAngle[][0]=Angles1[p]
			UnimodalMatrixAngle[][1]=Angles2[p]
		else
			angles1=angles2
			angles2=gnoise(sigma)
			UnimodalMatrixAngle[][count+1]=Angles2[p]
		endif
		ncontact=0
		for(n=0;n<9;n+=1)
			ncontact+=IsOverlapping(H,L,angles1[n],angles2[n])
		endfor
		count+=1
	while(ncontact>=minspokes)
	deletepoints/M=1 count,(400-count),UnimodalMatrixAngle
	print count
	return count
end

// FUNCTIONS FOR KYMOGRAPH FITTING (Figure S5)
function repeatedanalysis(basename,blockst,blocknum,trst,trbin,trnum,repnum)
string basename
variable blockst,blocknum,trst,trbin,trnum
variable repnum
variable i
make/O/N=(blocknum,repnum) maxsigmamatrix=0
for(i=0;i<repnum;i+=1)
VariableAverage(basename,blockst,blocknum)
wave mobilitymatrix,mobilitymatrixt

IsMatrixRandom(mobilitymatrixt,trst,trbin,trnum)
wave MaxSigmaDistance
maxsigmamatrix[][i]=MaxSigmaDistance[p]
endfor
matrixop AvMaxSigma=sumrows(maxsigmamatrix) 
AvMaxSigma/=repnum
duplicate/O maxsigmamatrix,MaxSigmaMatSssum
MaxSigmaMatSssum[][]=(maxsigmamatrix[p][q]-AvMaxSigma[p])^2
matrixop MaxSigmaStd=sumrows(MaxSigmaMatSssum) 
MaxSigmaStd/=repnum
End

function RandomBlock(peakposwave,blocksize)
wave peakposwave
variable blocksize
make/O/N=(dimsize(peakposwave,1)) CurrentBlockAv
make/N=(blocksize,dimsize(peakposwave,1))/FREE/O BlockPos,BlockDis
variable stmax,stx
stmax=dimsize(peakposwave,0)-blocksize-1
stx=(enoise(stmax)+stmax)/2 
BlockPos[][]=peakposwave[stx+p][q]
matrixop/FREE/O NewAvpos=sumcols(blockpos)^T
NewAvpos/=blocksize
BlockDis[][]=abs(BlockPos[p][q]-NewAvpos[q])
matrixop/O NewAvDisp=sumcols(BlockDis)^T
NewAvDisp/=blocksize
end

function VariableAverage(basename,blockst,blocknum)
string basename
variable blockst,blocknum
string matchname
matchname=basename+"*"
string listofwaves
listofwaves=WaveList(matchname, ";","" )
string currentwavename
variable j,i
make/O/N=(itemsinList(listofwaves),9,blocknum) MobilityMatrix
for(j=blockst;j<blockst+blocknum;j+=1)
	for(i=0;i<itemsinList(listofwaves);i+=1)
		currentwavename=stringfromlist(i,listofwaves)
		wave currentwave=$currentwavename
		RandomBlock(currentwave,j)
		wave NewAvDisp
		mobilitymatrix[i][][j-blockst]=NewAvDisp[q]
	endfor
endfor 
matrixop/o  mobilitymatrixt=mobilitymatrix^t
end

function IsMatrixRandom(mobilitymatrix,trst,trbin,trnum)
wave mobilitymatrix
variable trst,trbin,trnum
variable i
make/N=(dimsize(mobilitymatrix,0),dimsize(mobilitymatrix,1))/O CurrentMobMat
FractionBoundDifferentPerc(0,1/trnum,trnum,1000)
wave SimProb,SimFractionBound,SimEqualratios,SimStdRatios
//SimStdRatios/=sqrt(dimsize(mobilitymatrix,2))
make/O/N=(trnum)/FREE SigmaDistance
make/O/N=(dimsize(mobilitymatrix,2)) MaxSigmaDistance
for(i=0;i<dimsize(mobilitymatrix,2);i+=1)
	CurrentMobMat[][]=mobilitymatrix[p][q][i]
	SlidingTreshold(CurrentMobMat,trst,trbin,trnum)
	wave frbound,tres,equallyor
	make/O/N=(trnum) Xpoint=round(  (1-frbound[p])*trnum)
	SigmaDistance[]=abs((equallyor[p]-SimEqualratios[Xpoint[p]])/SimStdRatios[Xpoint[p]])
	MaxSigmaDistance[i]=wavemax(SigmaDistance)
endfor
end

function FindMaxInProfileNew(profilein,meanposwave,winsizewave)
wave profilein,meanposwave,winsizewave
variable i,j,stx,k
string newname,nameav,namedev,nameAvdev
newname=nameofwave(profilein)+"PeakPos"
nameav=nameofwave(profilein)+"AvPeakPos"
namedev=nameofwave(profilein)+"PeakDev"
nameAvdev=nameofwave(profilein)+"PeakAvDev"
make/O/N=9 $nameav,$nameAvdev
make/O/N=(dimsize(profilein,0),9) $newname,$namedev
wave avpos=$nameav
wave avdev=$nameAvdev
wave dev=$namedev
wave peak=$newname

for(j=1;j<=8;j+=1)//for each spoke
		stx=meanposwave[j]-winsizewave[j]/2
		make/O/FREE/N=(winsizewave[j]) Temp
	for(i=0;i<dimsize(profilein,0);i+=1) //for each timepoint
		temp[]=profilein[i][stx+p]
		wavestats temp
		peak[i][j-1]=stx+V_maxloc
	endfor	
endfor
make/O/FREE/N=(winsizewave[0]) Temp
		stx=meanposwave[0]-winsizewave[0]/2
for(i=0;i<dimsize(profilein,0);i+=1) //for each timepoint
		for(k=0;k<winsizewave[0];k+=1)
		if((stx+k)<dimsize(profilein,1))
			temp[k]=profilein[i][stx+k]
		else
			temp[k]=profilein[i][(stx+k)-dimsize(profilein,1)]
		endif
		endfor
		wavestats temp
		if((stx+V_maxloc)<dimsize(profilein,1))
			peak[i][8]=stx+V_maxloc
		else
			peak[i][8]=(stx+V_maxloc)-dimsize(profilein,1)
		endif
endfor	
make/O/N=(dimsize(profilein,0))/FREE temppeakpos,tempdev
for(j=0;j<8;j+=1)
	temppeakpos[]=peak[p][j]
	avpos[j]=mean(temppeakpos)
endfor

temppeakpos[]= (peak[p][j]<dimsize(profilein,1)/2) ? peak[p][j] : -(dimsize(profilein,1)-peak[p][j])
if(mean(temppeakpos)>0)
avpos[8]=mean(temppeakpos)
else
avpos[8]=dimsize(profilein,1)-mean(temppeakpos)
endif

dev[][0,7]=peak[p][q]-avpos[q]
if(avpos[8]<(dimsize(profilein,1)/2))
	dev[][8]= (peak[p][8]<dimsize(profilein,1)/2) ? peak[p][8]-avpos[8] : avpos[8]+(dimsize(profilein,1)-peak[p][8])
else
	dev[][8]= (peak[p][8]<dimsize(profilein,1)/2) ? dimsize(profilein,1)-avpos[8]+peak[p][8] : peak[p][8]-avpos[8]
endif


for(j=0;j<9;j+=1)
	tempdev[]=abs(dev[p][j])
	avdev[j]=mean(tempdev)
endfor
end

function FractionBoundDifferentPerc(probst,probbin,probnum,n)
variable probst,probbin,probnum,n
variable i,prob
make/O/N=(probnum) SimProb,SimFractionBound,SimEqualratios,SimStdRatios
for(i=0;i<probnum;i+=1)
prob=probst+i*probbin
RunMultipleSimulations(n,prob)
wave Simulatedratios,SimFrBound
SimProb[i]=prob
SimFractionBound[i]=mean(SimFrBound)
SimEqualratios[i]=mean(Simulatedratios)
SimStdRatios[i]= sqrt(variance(Simulatedratios))
endfor
end

Function SlidingTreshold(matrixin,trst,trbin,trnum)
wave matrixin
variable trst,trbin,trnum
variable i,tr
make/o/N=(trnum) frbound,tres,equallyor
make/FREE/N=(dimsize(matrixin,0),dimsize(matrixin,1)) Trmatrix
for(i=0;i<trnum;i+=1)
tr=trst+i*trbin
trmatrix[][] = matrixin[p][q]>tr ? 0 : 1
frbound[i]= sum(trmatrix)/numpnts(trmatrix)
tres[i]=tr
equallyor[i]=CountEquallyOriented(matrixin,tr)
endfor
end

function RunMultipleSimulations(numsim,prob)
variable numsim,prob
variable i
make/o/N=(numsim) Simulatedratios,SimFrBound
for(i=0;i<numsim;i+=1)
	SimulteAdherenceMatrix(20,prob)
	wave SimulatedAdMat
	Simulatedratios[i]=CountEquallyOriented(SimulatedAdMat,0.5)
	SimFrBound[i]=sum(SimulatedAdMat)/numpnts(SimulatedAdMat)
endfor
end

function CountEquallyOriented(matrixin,tres)
wave matrixin
variable tres
variable nighno=0,nn,equally=0,i,j,k
for(i=0;i<dimsize(matrixin,1);i+=1)//for each ring
	for(j=0;j<9;j+=1)
		for(k=0;k<2;k+=1)
			if(k==0)
				if((j-1)<0)
					nn=8
				else
					nn=j-1
				endif
			else
				if((j+1)>8)
					nn=0
				else
					nn=j+1
				endif
			endif
			nighno+=1
			if((matrixin[j][i]>tres) == (matrixin[nn][i]>tres))
				equally+=1
			endif
		endfor
	endfor
endfor
return equally/nighno
end

Function SimulteAdherenceMatrix(N,prob)
variable N,prob
make/N=(9,N)/o SimulatedAdMat
variable i,j
for(i=0;i<N;i+=1)
	for(j=0;j<9;j+=1)
		if((enoise(0.5)+0.5)>prob)
			SimulatedAdMat[j][i]=1
		else
			SimulatedAdMat[j][i]=0
		endif
	endfor
endfor
end



