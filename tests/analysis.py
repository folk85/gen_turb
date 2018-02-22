import numpy as np
import matplotlib.pyplot as plt

a = np.genfromtxt("c:/Users/sergey.m/work/2018/dns/35/user_io_stat.dat",skip_header=1)
fig, ax = plt.subplots()
j = -60
ax.plot(a[:j,0]*1000,a[:j,2]*1000, label=r'Mean Fl. Pos.')
ax.plot(a[:j,0]*1000,a[:j,3]*1000, label=r'Min Fl. Pos.')
ax.plot(a[:j,0]*1000,a[:j,4]*1000, label=r'Max Fl. Pos.')
ax.legend(loc=0)
# ax.set_xlim((0,500))
ax.set_xlabel(r'Time, ms')
ax.set_ylabel(r'Position, mm')
ax.grid()
plt.show()

x1 = np.linspace(a[0,0],a[j,0],70,endpoint=True)
b = np.interp(x1,a[:j,0],a[:j,2])
s = np.diff(b) / np.diff(x1)

fig, ax = plt.subplots()
j = -60
ax.plot(x1[1:]*1000,s, label=r'Mean Fl. Pos.')
# ax.plot(a[:j,0]*1000,a[:j,2]*1000, label=r'Mean Fl. Pos.')
# ax.plot(a[:j,0]*1000,a[:j,3]*1000, label=r'Min Fl. Pos.')
# ax.plot(a[:j,0]*1000,a[:j,4]*1000, label=r'Max Fl. Pos.')
ax.legend(loc=0)
# ax.set_xlim((0,500))
ax.set_xlabel(r'Time, ms')
ax.set_ylabel(r'Flame Velocity, m/s')
ax.grid()
plt.show()

class genvel():
    def __init__(self,name=None):
        if name !=None:
            coefs = np.genfromtxt(name,skip_header=1)
        else: 
            print("No files")
        self.acm = np.array([coefs[i,1:4] for i in range(len(coefs[:,0]))]) #.transpose
        self.asm = np.array([coefs[i,4:7] for i in range(len(coefs[:,0]))]) #.transpose
        self.bm = np.array([coefs[i,7:10] for i in range(len(coefs[:,0]))]) #.transpose
        self.cm = coefs[:,10]
        del(coefs)
        self.num = np.size(self.cm)
        print(np.shape(self.acm))
#         self.modes = np.dot()
#         self
    def get_vel(self,dx,dtime):
        db = np.matmul(self.bm,dx)
        db += np.dot(self.cm,dtime)
        dcos = np.cos(db)
        dsin = np.sin(db)
        res =  np.matmul(dcos,self.acm)
        res += np.matmul(dsin,self.asm)
        return res
#     def get_modes(self):
#         self.modes = 
class genvel2():
    def __init__(self,name=None):
        if name !=None:
            coefs = np.genfromtxt(name,skip_header=1)
        else: 
            print("No files")
        self.acm = coefs[:,1:4]
        self.asm = coefs[:,4:7]
        self.bm = np.transpose(coefs[:,7:10])
        self.cm = coefs[:,10]
        del(coefs)
        self.num = np.size(self.cm)
#         self.modes = np.dot()
#         self
        print(np.shape(self.acm),np.shape(self.bm))
    def get_vel(self,dx,dtime):
        db = np.matmul(dx,self.bm)
        db_add =  np.dot(self.cm,dtime)
        for i, _ in enumerate(db[:]):
            db[i,:] += db_add
        dcos = np.cos(db)
        dsin = np.sin(db)
        res =  np.matmul(dcos,self.acm)
        res += np.matmul(dsin,self.asm)
        return res
        
# print(coefs[0,1:])
# np.shape(acm)
# print(acm.shape,asm.shape,bm.shape,cm.shape)
# print(acm[0,:],asm[0,:],bm[0,:],cm[0])
# del(coefs)
# np.transpose
cl = genvel("c:/Users/sergey.m/work/2018/dns/35/user_coef.dat")
c2 = genvel2("c:/Users/sergey.m/work/2018/dns/35/user_coef.dat")


nx = [20]*3
dx = [5.0e-2/nx[0]]*3
dyy = [(np.arange(nx[i]) + .5 - nx[i]//2) * dx[i] for i in range(3)]
#dyy[1]
nel = np.size(dyy[0]) * np.size(dyy[1]) * np.size(dyy[2])
xp = np.zeros((nel,3))
# v = np.arange(50*3).reshape([3,50])
# d = np.arange(5*3).reshape([5,3])
# np.dot(d,v).shape
for i, ix in enumerate(dyy[0]):
    for j, jx in enumerate(dyy[1]):
        for k, kx in enumerate(dyy[2]):
            l = i * nx[1] * nx[2] + j * nx[2] + k
            xp[l,:] = np.array([dyy[0][i],dyy[1][j],dyy[2][k]])
            if (not(l%10000)): print(l)

vels = c2.get_vel(xp,1.0e-7)

print(vels[:,0].mean(),vels[:,1].mean(),vels[:,2].mean())
print(vels[:,0].std(),vels[:,1].std(),vels[:,2].std())
u_p = np.array([vels[:,i].std() for i in range(3)])
print(vels[:,:].std(),np.sqrt(np.sum(u_p**2)))