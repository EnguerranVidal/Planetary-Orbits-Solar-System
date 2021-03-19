# DATA MANIPULATION AND PLOTTING MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# FILE HANDLING MODULES
import os
import sys
import time


################################################################################
#--------- SIMULATION CLASS ---------#

class NBody_Engine():
    ''' This class serves as a calculation engine to solve the N-Body problem.'''
    def __init__(self):
        self.objects_name=[]
        self.objects_mass=[]
        self.objects_X=[]
        self.objects_V=[]
        self.objects_type=[]
        self.n_objects=0
        
    def define_objects(self,objects):
        ''' Defines and creates the arrays previously initiated using the inserted bodies' parameters '''
        n=len(objects)
        if self.n_objects==0:
            self.objects_X=np.zeros(shape=(n,3))
            self.objects_V=np.zeros(shape=(n,3))
            for i in range(n):
                self.objects_X[i]=objects[i][1]
                self.objects_V[i]=objects[i][2]
                self.objects_name.append(objects[i][0])
                self.objects_mass.append(objects[i][3])
        else:
            new_X=np.zeros(shape=(self.n_objects+n,3))
            new_V=np.zeros(shape=(self.n_objects+n,3))
            new_X[0:self.n_objects]=self.objects_X[0:self.n_objects]
            new_V[0:self.n_objects]=self.objects_V[0:self.n_objects]
            self.objects_X=new_X
            self.objects_V=new_V
            for i in range(n):
                self.objects_X[self.n_objects+i]=objects[i][1]
                self.objects_V[self.n_objects+i]=objects[i][2]
                self.objects_name.append(objects[i][0])
                self.objects_mass.append(objects[i][3])
                self.objects_type.append(objects[i][4])
        self.n_objects=self.n_objects+n     
 
    def objvect(self,objects_X,i,j):
        ''' Returns the vector from body i to j '''
        return self.objects_X[j]-self.objects_X[i]
    
    def objdist(self,objects_X,i,j):
        ''' Returns the distance between bodies i and j '''
        X=self.objvect(objects_X,i,j)
        return np.sqrt(X[0]**2+X[1]**2+X[2]**2)

    def gravconst(self):
        ''' the constant is in au^3/d^2/M_sol '''
        return 2.95912208286*10**(-4)
    
    def gravconst_SI(self):
        ''' the constant is in m^3/s^2/kg'''
        return 6.67408*10**(-11)
    
    def solar_mass(self):
        ''' the constant is in kg '''
        return 1.9884*10**(30)
    
    def astronomical_unit(self):
        ''' the constant is in meters '''
        return 1.49597870*10**(11)

    def acceleration(self,objects_X):
        ''' Returns an array containing the acceleration of the objects '''
        a=np.zeros_like(objects_X)
        n=self.n_objects
        for j in range(n):
            for k in range(n):
                if j!=k:
                    v=self.objvect(objects_X,j,k)
                    d=self.objdist(objects_X,j,k)
                    a[:][j]=a[:][j]+self.gravconst()*(self.objects_mass[k]*v)/(d**3)
        return a
    
    def compute(self,dt,method='Euler_explicit',focus_back=False):
        ''' Defines which computing method will be used according to the user's request '''
        if method=='Euler_explicit':
            self.compute_euler_explicit(dt)
        if method=='Euler_semi_implicit':
            self.compute_euler_semi_implicit(dt)
        if method=='Euler_symplectic':
            self.compute_euler_symplectic(dt)
        if method=='Heun':
            self.compute_Heun(dt)
        if method=='Runge_Kutta':
            self.compute_Runge_Kutta(dt)
        if focus_back==True:
            # Reference frame change : focus back on central body
            for i in range(self.n_objects):
                self.objects_X[i]=self.objects_X[i]-self.objects_X[0]
        
    def compute_euler_explicit(self,dt):
        ''' Calculates the system's next time step state using the explicit Euler method '''
        new_V=self.objects_V+dt*self.acceleration(self.objects_X)
        new_X=self.objects_X+dt*self.objects_V
        self.objects_X=new_X
        self.objects_V=new_V

    def compute_euler_semi_implicit(self,dt):
        ''' Calculates the system's next time step state using the semi-implicit Euler method '''
        new_V=self.objects_V+dt*self.acceleration(self.objects_X)
        new_X=self.objects_X+dt*new_V
        self.objects_X=new_X
        self.objects_V=new_V
    
    def compute_euler_symplectic(self,dt):
        ''' Calculates the system's next time step state using the symplectic Euler method '''
        new_X=self.objects_X+dt*self.objects_V
        new_V=self.objects_V+dt*self.acceleration(new_X)
        self.objects_X=new_X
        self.objects_V=new_V

    def compute_Heun(self,dt):
        ''' Calculates the system's next time step state using the Heun method '''
        k1_X=self.objects_V*dt
        k1_V=self.acceleration(self.objects_X)*dt
        k2_X=(self.objects_V+k1_V)*dt
        k2_V=self.acceleration(self.objects_X+k1_X)*dt
        new_X=self.objects_X+(k1_X+k2_X)/2
        new_V=self.objects_V+(k1_V+k2_V)/2
        self.objects_X=new_X
        self.objects_V=new_V

    def compute_Runge_Kutta(self,dt):
        ''' Calculates the system's next time step state using the Runge-Kutta method '''
        k1_X=self.objects_V*dt
        k1_V=self.acceleration(self.objects_X)*dt
        k2_X=(self.objects_V+k1_V/2)*dt
        k2_V=self.acceleration(self.objects_X+k1_X/2)*dt
        k3_X=(self.objects_V+k2_V/2)*dt
        k3_V=self.acceleration(self.objects_X+k2_X/2)*dt
        k4_X=(self.objects_V+k3_V)*dt
        k4_V=self.acceleration(self.objects_X+k3_X)*dt
        new_X=self.objects_X+(k1_X+2*k2_X+2*k3_X+k4_X)/6
        new_V=self.objects_V+(k1_V+2*k2_V+2*k3_V+k4_V)/6
        self.objects_X=new_X
        self.objects_V=new_V

################################################################################
#--------- PLANETARY SYSTEM CLASS ---------#

class Planetary_System():
    ''' Main class directing the entire program.'''
    def __init__(self,ephemeride_file):
        self.current_dir=os.path.dirname(os.path.abspath(__file__))
        self.ephemeride_file=ephemeride_file
        self.engine=NBody_Engine()
        self.database=Ephemeride_Database(self.ephemeride_file)
        objects=self.database.load_data()
        self.engine.define_objects(objects)
        self.time=0
        self.is_new=True
        self.saves_file=None
        self.n_saves=0

    def new_session(self,ephemeride_file):
        ''' Initializes a new session, in order to compute new values, starting from scratch '''
        self.engine=NBody_Engine()
        self.ephemeride_file=ephemeride_file
        self.database=Ephemeride_Database(self.ephemeride_file)
        objects=self.database.load_data()
        self.engine.define_objects(objects)
        self.time=0
        self.is_new=True
        self.saves_file=None
        self.n_saves=0

    def load_session(self,save_file_name):
        ''' Loads a previous session, in order to use older, already computed data '''
        new_path=self.current_dir+"\\logs\\"+save_file_name
        assert os.path.exists(new_path)==True,"ERROR : File does not exists, try a different name."
        self.saves_file=save_file_name
        self.is_new=False
        [x,m]=self.load_save_info()
        if self.ephemeride_file!=x:
            print("WARNING : ephemeride file initialized does not match with the one in save file.")
            print("Replacing "+self.ephemeride_file+'data by '+x+' data .....')
            self.ephemeride_file=x
            self.database=Ephemeride_Database(self.ephemeride_file)
            objects=self.database.load_data()
            self.engine.define_objects(objects)
        with open(new_path,"r") as file:
            file.readline()
            data_needed=[0,1,2,3,4,5,6]
            for i in range(m):
                data=self.load_state(file,data_needed)
            X=np.zeros_like(self.engine.objects_X)
            V=np.zeros_like(self.engine.objects_V)
            X[:,0]=data[1]
            X[:,1]=data[2]
            X[:,2]=data[3]
            V[:,0]=data[4]
            V[:,1]=data[5]
            V[:,2]=data[6]
            self.engine.objects_X=X
            self.engine.objects_V=V
            self.time=data[0]
        self.n_saves=m
        
    def save_state(self):
        ''' Saves the computed data in the file containing data from older time steps '''
        self.n_saves=self.n_saves+1
        parameters=np.zeros(shape=(self.engine.n_objects,6))
        for i in range(1,self.engine.n_objects):
            parameters[i]=orbital_parameters(self.engine.objects_X[i],self.engine.objects_V[i],
                                             self.engine.gravconst(),self.engine.objects_mass[0])
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        X,Y,Z,Xp,Yp,Zp,a,e,i,Omega,w,theta=[],[],[],[],[],[],[],[],[],[],[],[]
        for j in range(self.engine.n_objects):
            X.append(str(self.engine.objects_X[j,0]))
            Y.append(str(self.engine.objects_X[j,1]))
            Z.append(str(self.engine.objects_X[j,2]))
            Xp.append(str(self.engine.objects_V[j,0]))
            Yp.append(str(self.engine.objects_V[j,1]))
            Zp.append(str(self.engine.objects_V[j,2]))
            a.append(str(parameters[j,0]))
            e.append(str(parameters[j,1]))
            i.append(str(parameters[j,2]))
            Omega.append(str(parameters[j,3]))
            w.append(str(parameters[j,4]))
            theta.append(str(parameters[j,5]))           
        with open(new_path,'a') as file:
            file.write(str(self.time)+'\n')
            file.write(' '.join(X)+'\n')
            file.write(' '.join(Y)+'\n')
            file.write(' '.join(Z)+'\n')
            file.write(' '.join(Xp)+'\n')
            file.write(' '.join(Yp)+'\n')
            file.write(' '.join(Zp)+'\n')
            file.write(' '.join(a)+'\n')
            file.write(' '.join(e)+'\n')
            file.write(' '.join(i)+'\n')
            file.write(' '.join(Omega)+'\n')
            file.write(' '.join(w)+'\n')
            file.write(' '.join(theta)+'\n')
            file.close()
            
    def load_state(self,file,data_needed):
        ''' Loads the computed data for a given time step from the open session file '''
        assert self.is_new==False,"ERROR : The save file has yet to be rendered, try doing calculations first."
        data=[]
        holder=file.readline()
        if 0 in data_needed:
            data.append(float(holder))
        for i in range(1,13):
            holder=file.readline()
            if i in data_needed:
                data.append(str_to_float_list(holder))
        return data

    def load_save_info(self):
        ''' Returns information on the specifics of which file has been loaded '''
        assert self.is_new==False,"ERROR : The save file has yet to be rendered, try doing calculations first."
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        with open(new_path,'r') as file:
            initial_line=file.readline()
            file.close()
        initial_line=initial_line.split()
        initial_line.pop(0)
        initial_line[1]=int(initial_line[1])
        return initial_line
            
    def RUN(self,dt,T,skip,method='Euler_explicit'):
        ''' Runs the calculations, using the ephemeride data initialized, and saves them at each time step '''
        if self.is_new==True:
            logs_path=self.current_dir+"\\logs"
            if os.path.exists(logs_path)==False:
                os.mkdir(logs_path)
            self.saves_file=session_name()
            new_path=self.current_dir+"\\logs\\"+self.saves_file
            with open(new_path,'w') as file:
                file.write('Base_File '+self.ephemeride_file+' '+str(self.n_saves)+'\n')
                file.close()
        self.save_state()
        initial_time=self.time
        last_snap_time=self.time
        print("Beginning Calculations")
        while self.time<initial_time+T:
            self.engine.compute(dt,method=method,focus_back=True)
            self.time=self.time+dt
            if self.time-last_snap_time>=skip:
                self.save_state()
                last_snap_time=self.time
        print("Calculations Finished")
        # Updating the number of snapshots contained inside the file
        file=open(self.current_dir+"\\logs\\"+self.saves_file,"r")
        lines=file.readlines()
        lines[0]='Base_File '+self.ephemeride_file+' '+str(self.n_saves)+'\n'
        file.close()
        file=open(self.current_dir+"\\logs\\"+self.saves_file,"w")
        file.writelines(lines)
        file.close()
        self.is_new=False

    def display_3D(self,labels=True):
        ''' Displays a 3D animation of the planetary system, with the star at the center and the planets' orbit around it using the computed data '''
        assert self.is_new==False,"ERROR : Cannot display System since no calculations have taken place."
        print("Displaying 3D System")
        # Creating the 3D figure
        fig=plt.figure(figsize=(12,12))
        ax=fig.gca(projection='3d')
        ax.set_title('t = 0.0 days')
        ax.set_xlim3d(-50,50)
        ax.set_ylim3d(-50,50)
        ax.set_zlim3d(-50,50)
        xLabel=ax.set_xlabel('\nX [ au ]',linespacing=3.2)
        yLabel=ax.set_ylabel('\nY [ au ]',linespacing=3.1)
        zLabel=ax.set_zlabel('\nZ [ au ]',linespacing=3.4)
        # Accessing the save file
        [x,m]=self.load_save_info()
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        initial_needed_data=[0,1,2,3,4,5,6,7,8,9,10,11,12]
        with open(new_path,'r') as file:
            file.readline() #Skipping first line
            data=self.load_state(file,initial_needed_data)
            graph=ax.scatter(data[1],data[2],data[3],c='y',edgecolor="k")
            orbits=[]
            for i in range(1,self.engine.n_objects):
                X,Y,Z=find_trajectory(data[7][i],data[8][i],data[9][i],
                                      data[10][i],data[11][i],data[12][i],120)
                orbits.append(ax.plot(X,Y,Z))
            if labels==True:
                Labels=[]
                for i in range(0,self.engine.n_objects):
                    Labels.append(ax.text(data[1][i],data[2][i],data[3][i],self.engine.objects_name[i], (1,1,1)))
            fig.show()
            plt.pause(3)
            for j in range(1,m):
                plt.pause(0.04)
                data=self.load_state(file,initial_needed_data)
                graph._offsets3d=(data[1],data[2],data[3])
                for k in range(1,self.engine.n_objects):
                    X,Y,Z=find_trajectory(data[7][k],data[8][k],data[9][k],
                                      data[10][k],data[11][k],data[12][k],120)
                    line=orbits[k-1][0]
                    line.set_data(X,Y)
                    line.set_3d_properties(Z)
                    orbits[k-1]=[line]
                ax.set_title('t = '+str(round(data[0],))+' days')
                if labels==True:
                    for i in range(0,self.engine.n_objects):
                        Labels[i].set_position((data[1][i],data[2][i]))
                        Labels[i].set_3d_properties(data[3][i],(1,1,1))
                plt.draw()
            file.close()
    
    def apsidal_precession(self,displayed='all'):
        ''' Calculates and displays the apsidal precession of the bodies requested by the user over time '''
        assert self.is_new==False,"ERROR : The save file has yet to be rendered, try doing calculations first."
        if displayed=='all':
            planets=[]
            for i in range(1,self.engine.n_objects):
                planets.append(i)
        else:
            if type(displayed)==type('n'):
                assert displayed in self.engine.objects_name,displayed+" is not in the Ephemeride file objects list."
                planets=[self.engine.objects_name.index(displayed)]
            if type(displayed)==type(['n']):
                planets=[]
                for name in displayed:
                    assert name in self.engine.objects_name,name+" is not in the Ephemeride file objects list."
                    planets.append(self.engine.objects_name.index(name))
        if len(planets)>1:
            print(" Displaying apsidal precessions")
        else:
            print(" Displaying apsidal precession")
        # Accessing the save file
        [x,m]=self.load_save_info()
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        initial_needed_data=[0,11]
        n=self.engine.n_objects
        Time=np.zeros(shape=(m,1))
        w0=np.zeros(shape=(1,n))
        w=np.zeros(shape=(m,n))
        with open(new_path,'r') as file:
            file.readline() #Skipping first line
            data=self.load_state(file,initial_needed_data)
            w0=data[1]
            precession=np.zeros(shape=(m,n))
            Time[0]=data[0]
            for i in range(1,m):
                data=self.load_state(file,initial_needed_data)
                w[i]=data[1]
                precession[i]=w[i]-w0
                Time[i]=data[0]
                precession[i]=precession[i]*180/np.pi #from rad to deg
        file.close()
        for j in planets:
            plt.plot(Time,precession[:,j], label=self.engine.objects_name[j])
        plt.title('Apsidal precession over time')
        plt.xlabel('Time (Days)')
        plt.ylabel('Orbital Shift (°)')
        plt.legend()    
        plt.show()
          
    
    def display_perihelion(self,displayed='all'):
        ''' Calculates and displays the perihelion of the bodies requested by the user over time '''
        assert self.is_new==False,"ERROR : The save file has yet to be rendered, try doing calculations first."
        if displayed=='all':
            planets=[]
            for i in range(1,self.engine.n_objects):
                planets.append(i)
        else:
            if type(displayed)==type('n'):
                assert displayed in self.engine.objects_name,displayed+" is not in the Ephemeride file objects list."
                planets=[self.engine.objects_name.index(displayed)]
            if type(displayed)==type(['n']):
                planets=[]
                for name in displayed:
                    assert name in self.engine.objects_name,name+" is not in the Ephemeride file objects list."
                    planets.append(self.engine.objects_name.index(name))
        if len(planets)>1:
            print(" Displaying perihelions")
            plt.title('Perihelions values')
        else:
            print(" Displaying perihelion")
            plt.title('Perihelion values')
        print(planets)
        # Accessing the save file
        [x,m]=self.load_save_info()
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        initial_needed_data=[0,7,8]
        n=self.engine.n_objects
        Times=np.zeros(shape=(m,))
        A=np.zeros(shape=(m,n))
        E=np.zeros(shape=(m,n))
        with open(new_path,'r') as file:
            file.readline() #Skipping first line
            for i in range(m):
                data=self.load_state(file,initial_needed_data)
                Times[i]=data[0]
                A[i]=data[1]
                E[i]=data[2]
            file.close()
        R=A*(1-E)
        for j in planets:
            plt.plot(Times,R[:,j],label=self.engine.objects_name[j])
        plt.legend()
        plt.xlabel("Time (Days)")
        plt.ylabel("Periapsis (AUs)")
        plt.show()

    def energy_conservation(self):
        '''Calculates the mechanical energy of the entire system at each time and plots it over time.'''
        assert self.is_new==False,"ERROR : The save file has yet to be rendered, try doing calculations first."
        [x,m]=self.load_save_info()
        new_path=self.current_dir+"\\logs\\"+self.saves_file
        initial_needed_data=[0,7]
        n=self.engine.n_objects
        Times=np.zeros(shape=(m,))
        E=np.zeros(shape=(m,n))
        with open(new_path,'r') as file:
            file.readline() #Skipping first line
            for i in range(m):
                data=self.load_state(file,initial_needed_data)
                Times[i]=data[0]
                energy=0
                for j in range(1,n):
                    mj=self.engine.objects_mass[j]
                    M=self.engine.objects_mass[0]
                    G=self.engine.gravconst()
                    energy=energy-mj*M*G/(2*data[1][j])
                E[i]=energy
            file.close()
        plt.plot(Times,E)
        plt.xlabel("Time (Days)")
        plt.ylabel("Energy")
        plt.title("Total System Energy")
        plt.show()
        


################################################################################
#--------- EPHEMERIDE CLASS ---------#

class Ephemeride_Database():
    ''' Initializes and handles the ephemeride data on the planetary system's initial conditions '''
    def __init__(self,ephemeride_file):
        self.current_dir=os.path.dirname(os.path.abspath(__file__))
        self.ephemeride_file=ephemeride_file
        self.path=self.current_dir+"\\ephem\\"
        self.filename=self.current_dir+"\\ephem\\"+ephemeride_file
        assert os.path.exists(self.path)==True,"The ephemerides folder is not present, please create it using the name 'ephem'."
        assert os.path.exists(self.filename)==True,self.ephemeride_file+" is not in the ephemerides folder."
        with open(self.filename,'r') as file:
            objects=file.readlines()
        self.reference_time=objects[0]
        self.labels=objects[1]
        objects.pop(0)
        objects.pop(0)
        n=len(objects)
        self.catalogue=objects

    def add_object(self):
        ''' Adds an object and its initial parameters in the file containing such data on other objects '''
        name=str(input("Object's name :"))
        mass=float(input("Object's mass :"))
        x=float(input("Object's x :"))
        y=float(input("Object's y :"))
        z=float(input("Object's z :"))
        xp=float(input("Object's xp :"))
        yp=float(input("Object's yp :"))
        zp=float(input("Object's zp :"))
        string='\n'+name+','+str(mass)+','+str(x)+','+str(y)+','+str(z)+','
        string=string+str(xp)+','+str(yp)+','+str(zp)
        with open(self.filename,'a') as file:
            file.write(string)
            file.close()

    def load_data(self):
        ''' Loads the ephemeride data from the given file '''
        objects=[]
        n=len(self.catalogue)
        for i in range(n):
            object_i=self.catalogue[i].split(',')
            m=len(object_i)
            name=object_i[0]
            mass=float(object_i[1])
            x=float(object_i[2])
            y=float(object_i[3])
            z=float(object_i[4])
            xp=float(object_i[5])
            yp=float(object_i[6])
            zp=float(object_i[7])
            object_i=[name,[x,y,z],[xp,yp,zp],mass]
            objects.append(object_i)
        return objects

 
################################################################################
#--------- FUNCTIONS ---------#

def quadrant(cos_i,sin_i):
    ''' Returns the real angle by using the four quadrants in trigonometry '''
    if cos_i>=0 and sin_i>=0:  # Quadrant 1
        return np.arccos(cos_i)
    if cos_i<0 and sin_i>=0:  # Quadrant 2
        return np.pi-np.arccos(np.abs(cos_i))
    if cos_i>=0 and sin_i<0:  # Quadrant 4
        return 2*np.pi-np.arccos(cos_i)
    if cos_i<0 and sin_i<0:  # Quadrant 3
        return np.pi+np.arccos(np.abs(cos_i))
        

def orbital_parameters(R,V,G,M):
    ''' Calculates the orbital parameters used to describe a Keplerian orbit '''
    r=np.linalg.norm(R)
    v=np.linalg.norm(V)
    energy=(v**2)/2-(G*M)/r
    a=-(G*M)/(2*energy)
    E=(1/(G*M))*((v**2-(G*M)/r)*R-np.dot(R,V)*V)
    e=np.linalg.norm(E)
    H=np.cross(R,V)
    h=np.linalg.norm(H)
    K=np.array([0,0,1])
    I=np.array([1,0,0])
    J=np.array([0,1,0])
    i=np.arccos(np.dot(K,H)/h)
    N=np.cross(K,H)
    n=np.linalg.norm(N)
    cos_omega=np.dot(I,N)/n
    sin_omega=np.dot(J,N)/n
    omega=quadrant(cos_omega,sin_omega)
    if E[2]>=0:
        w=np.arccos((np.dot(N,E))/(n*e))
    else:
        w=2*np.pi-np.arccos((np.dot(N,E))/(n*e))
    if np.dot(R,V)>=0:
        theta=np.arccos((np.dot(E,R))/(e*r))
    else:
        theta=2*np.pi-np.arccos((np.dot(E,R))/(e*r))
    return [a,e,i,omega,w,theta]


def find_trajectory(a,e,i,Omega,w,theta,N):
    ''' Computes the trajectory of a body by using its Keplerian orbital parameters '''
    rot_Omega=np.array([[np.cos(Omega),np.sin(Omega),0],
                        [-np.sin(Omega),np.cos(Omega),0],
                        [0,0,1]])
    rot_w=np.array([[np.cos(w),np.sin(w),0],
                    [-np.sin(w),np.cos(w),0],
                    [0,0,1]])
    rot_i=np.array([[1,0,0],
                    [0,np.cos(i),np.sin(i)],
                    [0,-np.sin(i),np.cos(i)]])
    rotation_matrix=np.matmul(rot_w,np.matmul(rot_i,rot_Omega))
    if e<1: #Ellipse case
        Thetas=np.linspace(-np.pi,np.pi,num=N)
        Radiuses=(a*(1-e**2))/(1+e*np.cos(Thetas))
        Xs=[]
        Ys=[]
        Zs=[]
        n=len(Thetas)
        for i in range(n):
            X=np.array([Radiuses[i]*np.cos(Thetas[i]),Radiuses[i]*np.sin(Thetas[i]),0])
            X=np.matmul(np.linalg.inv(rotation_matrix),X)
            Xs.append(X[0])
            Ys.append(X[1])
            Zs.append(X[2])
        return [np.array(Xs),np.array(Ys),np.array(Zs)]
    if e>=1: #Parabola/hyperbola case
        limit=np.arccos(-1/e)
        Thetas=np.linspace(-limit,limit,num=N)
        Radiuses=(a*(1-e**2))/(1+e*np.cos(Thetas))
        Xs=[]
        Ys=[]
        Zs=[]
        n=len(Thetas)
        for i in range(n):
            X=np.array([Radiuses[i]*np.cos(Thetas[i]),Radiuses[i]*np.sin(Thetas[i]),0])
            X=np.matmul(np.linalg.inv(rotation_matrix),X)
            Xs.append(X[0])
            Ys.append(X[1])
            Zs.append(X[2])
        return [np.array(Xs),np.array(Ys),np.array(Zs)]

def str_to_float_list(string):
    L=string.split(' ')
    n=len(L)
    for i in range(n):
        L[i]=float(L[i])
    return L

def kepler_to_cartesian(a,e,i,Omega,w,theta,star_mass,planet_mass):
    ''' Gives the position and speed vectors in a stellarcentric reference frame.
     a needs to be in AUs the main angles in radians, star_mass in units 
    of solar mass and planet_mass in units of Jupiter's mass. '''
    jupiter_mass=1.898*10**27 #kg
    sun_mass=1.989*10**30 #kg
    G=2.95912208286*10**(-4) #ua^3/d^2/M_sol
    p=a*(1-e**2)
    r=p/(1+e*np.cos(theta))
    planet_mass=planet_mass*jupiter_mass/sun_mass
    mu=G*star_mass
    h=np.sqrt(p*mu)
    R=np.array([r*np.cos(theta),r*np.sin(theta),0])
    V=np.array([-mu*np.sin(theta)/h,mu*(e+np.cos(theta))/h,0])
    # Rotation Matrix
    rot_Omega=np.array([[np.cos(Omega),np.sin(Omega),0],
                        [-np.sin(Omega),np.cos(Omega),0],
                        [0,0,1]])
    rot_w=np.array([[np.cos(w),np.sin(w),0],
                    [-np.sin(w),np.cos(w),0],
                    [0,0,1]])
    rot_i=np.array([[1,0,0],
                    [0,np.cos(i),np.sin(i)],
                    [0,-np.sin(i),np.cos(i)]])
    rotation_matrix=np.matmul(rot_w,np.matmul(rot_i,rot_Omega))
    # Transformation of Position and Speed
    R=np.matmul(np.linalg.inv(rotation_matrix),R)
    V=np.matmul(np.linalg.inv(rotation_matrix),V)
    print("Object's mass (in M_sol) : ",planet_mass)
    print("Object's X,Y,Z (in M_sol) : "+str(R[0])+','+str(R[1])+','+str(R[2]))
    print("Object's Xp,Yp,Zp (in M_sol) : "+str(V[0])+','+str(V[1])+','+str(V[2]))

def session_name():
    ''' Names the files with the given date and time format: yyyy-mm-dd-hours-mins-secs '''
    t0=time.time()
    struct=time.localtime(t0)
    string=str(struct.tm_year)+'-'
    # MONTHS
    n_months=str(struct.tm_mon)
    if len(n_months)==1:
        n_months='0'+n_months
    string=string+n_months+'-'
    # DAYS
    n_days=str(struct.tm_mday)
    if len(n_months)==1:
        n_days='0'+n_days
    string=string+n_days+'-'
    # HOURS
    n_hours=str(struct.tm_hour)
    if len(n_hours)==1:
        n_hours='0'+n_hours
    string=string+n_hours+'-'
    # MINUTES
    n_mins=str(struct.tm_min)
    if len(n_mins)==1:
        n_mins='0'+n_mins
    string=string+n_mins+'-'
    # SECONDS
    n_secs=str(struct.tm_sec)
    if len(n_secs)==1:
        n_secs='0'+n_secs
    string=string+n_secs+'.txt'
    return string
