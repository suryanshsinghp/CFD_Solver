
import numpy as np
import matplotlib.pyplot as plt
filename='drag_lift.dat'
filenameppe='PPE_residual.dat'
filenameade='ADE_residual.dat'
fileprobep='pressprobe.dat'
fileprobeuv='velprobe.dat'
filesurfpress='surfacepress.dat'
steadyt=15
iter_col=1
res_col=2
time_col=1
drag_col=2
lift_col=3
skip_ld=10
uprobe_col=2
vprobe_col=3
pprobe_col=2
tprobe_col=1

iter_ppe= np.genfromtxt(filenameppe,skip_header=0,usecols=iter_col-1)
iter_ade= np.genfromtxt(filenameade,skip_header=0,usecols=iter_col-1)
res_ppe= np.genfromtxt(filenameppe,skip_header=0,usecols=res_col-1)
res_ade= np.genfromtxt(filenameade,skip_header=0,usecols=res_col-1)
time=np.genfromtxt(filename,skip_header=skip_ld,usecols=time_col-1)
drag=np.genfromtxt(filename,skip_header=skip_ld,usecols=drag_col-1)
lift=np.genfromtxt(filename,skip_header=skip_ld,usecols=lift_col-1)
uprobe=np.genfromtxt(fileprobeuv,skip_header=skip_ld,usecols=uprobe_col-1)
vprobe=np.genfromtxt(fileprobeuv,skip_header=skip_ld,usecols=vprobe_col-1)
pprobe=np.genfromtxt(fileprobep,skip_header=skip_ld,usecols=pprobe_col-1)
tprobe=np.genfromtxt(fileprobep,skip_header=skip_ld,usecols=tprobe_col-1)



plt.plot(iter_ppe,res_ppe, color='dodgerblue',label="Residuals PPE")
plt.plot(iter_ade,res_ade, color='k',label="Residuals ADE")
plt.title("Global Residual Tracking")
plt.legend(bbox_to_anchor=(0.50, 0.50))
plt.xlabel(r'Iteration',fontsize=15)
plt.ylabel(r'Residual',fontsize=15)


file='residual.eps'
plt.savefig(file)
plt.close()

plt.plot(time,drag, color='dodgerblue',label="$c_D$")
plt.plot(time,lift, color='k',label="$c_L$")
plt.title("Lift and Drag coefficient")
plt.legend()
plt.xlabel(r'Time',fontsize=15)
plt.ylabel(r'Lift/Drag coefficient',fontsize=15)

file='draglift.eps'
plt.savefig(file)
plt.close()

plt.plot(tprobe,uprobe, color='dodgerblue',label="u")
plt.plot(tprobe,vprobe, color='k',label="v")
plt.plot(tprobe,pprobe,color='firebrick',label="pressure")
plt.title("Velocity and Pressure at probe location")
plt.legend()
plt.xlabel(r'Time',fontsize=15)
plt.ylabel(r'velocity/pressure',fontsize=15)

file='probe.eps'
plt.savefig(file)
plt.close()

angle=np.genfromtxt(filesurfpress,skip_header=0,usecols=1-1)
surfp=np.genfromtxt(filesurfpress,skip_header=0,usecols=2-1)
arr1inds = angle.argsort()
angle = angle[arr1inds[::-1]]
surfp = surfp[arr1inds[::-1]]
plt.plot(angle,surfp, color='k')      
plt.title("Cylinder Surface Pressure")
plt.xlabel(r'Angle',fontsize=15)
plt.ylabel(r'Pressure',fontsize=15)
file='surfpres.eps'
plt.savefig(file)
plt.close()


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
ststep=find_nearest(time,steadyt)
meandrag=np.mean(drag[ststep:])
meanlift=np.mean(lift[ststep:])
print("mean drag: ",meandrag)
print("mean lift: ",meanlift)


spectrum=np.fft.fft(vprobe[ststep:np.shape(tprobe)[0]])
freq = np.fft.fftfreq(len(spectrum))
freq=freq/0.001
#plt.plot(freq, abs(spectrum))
threshold = 0.5 * max(abs(spectrum))
mask = abs(spectrum) > threshold
peaks = freq[mask]
print("shedding frequency: ",max(peaks))

plt.plot(freq[0:80], abs(spectrum[0:80]))     
plt.title("FFT of v-velocity probe")
plt.xlabel(r'frequency',fontsize=15)
file='shedfreq.eps'
plt.savefig(file)
plt.close()