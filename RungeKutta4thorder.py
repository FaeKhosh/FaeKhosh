
import sys
import numpy as np
import matplotlib.pyplot as plt   

def dY1(y2):
    return y2
 
def dY2(y3):
    return y3
 
def dY3(y1,y3,y4,y5,cMu,Tinf):
    return - y3*((y5/(2*y4)) - (y5/(y4+(cMu/Tinf))) ) - (y1*y3)*( (y4+(cMu/Tinf)) / (np.sqrt(y4)*(1+(cMu/Tinf)) ))
    
def dY4(y5):
    return y5
 
def dY5(y1, y3, y4, y5, cMu, Tinf, Minf,Pr,gamma):
    return -(y5**2)*((1/(2*y4)) - (1/(y4+(cMu/Tinf))) ) - Pr *(y1*y5/np.sqrt(y4))*( (y4+(cMu/Tinf))/(1+(cMu/Tinf)) ) -(gamma - 1)* Pr*(Minf**2)*(y3**2)
# =============================================================================
def RK(N, h,y1,y2,y3,y4,y5,cMu,Tinf,Pr,gamma, Minf):
    for i in range(N):
        k1_y1 = dY1(y2[i]) 
        k1_y2 = dY2(y3[i])
        k1_y3 = dY3(y1[i], y3[i], y4[i], y5[i],cMu,Tinf)  
        k1_y4 = dY4(y5[i]) 
        k1_y5 = dY5(y1[i], y3[i], y4[i], y5[i],cMu,Tinf, Minf,Pr,gamma) 
        
        k2_y1 = dY1(y2[i]+ (h*0.5*k1_y2))
        k2_y2 = dY2(y3[i]+ (h*0.5*k1_y3))
        k2_y3 = dY3(y1[i]+(h*0.5*k1_y1) , y3[i]+(h*0.5*k1_y3), y4[i]+(h*0.5*k1_y4), y5[i]+(h*0.5*k1_y5),cMu,Tinf) 
        k2_y4 = dY4(y5[i]+ (h*0.5*k1_y5))
        k2_y5 = dY5(y1[i]+(h*0.5*k1_y1) , y3[i]+(h*0.5*k1_y3), y4[i]+(h*0.5*k1_y4), y5[i]+(h*0.5*k1_y5),cMu,Tinf, Minf,Pr,gamma) 
               
        k3_y1 = dY1(y2[i]+ (h*0.5*k2_y2))
        k3_y2 = dY2(y3[i]+ (h*0.5*k2_y3))
        k3_y3 = dY3(y1[i]+(h*0.5*k2_y1), y3[i]+(h*0.5*k2_y3), y4[i]+(h*0.5*k2_y4), y5[i]+(h*0.5*k2_y5),cMu,Tinf) 
        k3_y4 = dY4(y5[i]+ (h*0.5*k2_y5))
        k3_y5 = dY5(y1[i]+(h*0.5*k2_y1), y3[i]+(h*0.5*k2_y3), y4[i]+(h*0.5*k2_y4), y5[i]+(h*0.5*k2_y5),cMu,Tinf, Minf,Pr,gamma) 
        
        k4_y1 = dY1(y2[i]+ (h*0.5*k3_y2))
        k4_y2 = dY2(y3[i]+ (h*0.5*k3_y3))
        k4_y3 = dY3(y1[i]+(h*k3_y1), y3[i]+(h*k3_y3), y4[i]+(h*k3_y4), y5[i]+(h*k3_y5),cMu,Tinf) 
        k4_y4 = dY4(y5[i]+ (h*0.5*k3_y5))
        k4_y5 = dY5(y1[i]+(h*k3_y1), y3[i]+(h*k3_y3), y4[i]+(h*k3_y4), y5[i]+(h*k3_y5),cMu,Tinf, Minf,Pr,gamma) 
        
        y1.append(y1[i] + ((h/6)* (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1 ) ) )
        y2.append(y2[i] + ((h/6)* (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2 ) ) )
        y3.append(y3[i] + ((h/6)* (k1_y3 + (2*k2_y3) + (2*k3_y3) + k4_y3))  )
        y4.append(y4[i] + ((h/6)* (k1_y4 + 2*k2_y4 + 2*k3_y4 + k4_y4 ) ) )
        y5.append(y5[i] + ((h/6)* (k1_y5 + (2*k2_y5) + (2*k3_y5) + k4_y5)) )

    return y1,y2,y3,y4,y5

def selfsimilar(Minf, Tinf=273.15, Tw=1.007, eta_max=15, N=1000, eps=1e-9,M= 3):    
    delta_eta= eta_max/N
    gamma = 1.4  
    cMu = 110.4  
    Pr = 0.72
    delta = 1e-2
    eta  = [i*delta_eta for i in range(N)]
    adi == 1 #no specific boundary condition 
    
    # Initializing the solution vectors:
    y1 = [0] # f,   
    y2 = [0] # f'
    y3 = [0] # f''
    y4 = [0] # ρ(η)
    y5 = [0] # ρ(η)'
    if adi == 1:
        y5[0] = 0 
        alpha = 0.4700  # Initial Guess
        beta = 1.12      # Initial Guess

    y2n = 1 # f(infinity) = 1
    y4n = 1 # rho(infinity) = 1
    for indx in range(M):
        if adi==1:
            y1[0] = 0
            y2[0] = 0
            y5[0] = 0 
            y3[0] = alpha # Initial Guess
            y4[0] = beta # Initial Guess 
		# Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,delta_eta,y1,y2,y3,y4,y5,cMu,Tinf,Pr,gamma,Minf)

        y2_zero = y2[-1]  
        y4_zero = y4[-1]  
        
        # Initializing the solution vectors
        y1 = [0]   # f
        y2 = [0]   # f'
        y3 = [0]   # f''
        y4 = [0]   # ρ(η)
        y5 = [0]   # ρ(η)'

        # Small number addition for Newton's iteration method
        if adi==1:
            y1[0] = 0 
            y2[0] = 0 
            y5[0] = 0 
            y3[0] = alpha + delta   # Initial Guess + Small number
            y4[0] = beta   # Initial Guess

		# Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,delta_eta,y1,y2,y3,y4,y5,cMu,Tinf,Pr,gamma,Minf)
        # Storing the freestream values for Newton's iteration method
        y2n1 = y2[-1] 
        y4n1 = y4[-1]  
        # Small number addition for Newton's iteration method
        y1 = [0]   # f
        y2 = [0]   # f'
        y3 = [0]   # f''
        y4 = [0]   # ρ(η)
        y5 = [0]   # ρ(η)'
        if adi==1:
            y1[0] = 0
            y2[0] = 0 
            y5[0] = 0 
            y3[0] = alpha   # Initial Guess
            y4[0] = beta + delta  # Initial Guess + Small number            

		# Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,delta_eta,y1,y2,y3,y4,y5,cMu,Tinf,Pr,gamma, Minf)
        # Storing the freestream values for Newton's iteration method
        y2n2 = y2[-1]  
        y4n2 = y4[-1]  

        # Calculation of the next initial guess with Newton's iteration method
        p11 = (y2n1-y2_zero)/delta
        p21 = (y4n1-y4_zero)/delta
        p12 = (y2n2-y2_zero)/delta
        p22 = (y4n2-y4_zero)/delta
        
        r1 = 1-y2_zero
        r2 = 1-y4_zero

        A = np.array([[p11, p12], [p21, p22]])
        B = np.array([r1, r2])
        delta_alpha= np.linalg.inv(A).dot(B)[0]
        delta_beta =  np.linalg.inv(A).dot(B)[1]

        alpha = alpha + delta_alpha
        beta = beta + delta_beta

 
        if (abs(y2[-1]-1.0)< eps) and (abs(y4[-1]-1.0)<eps) and (abs(y2n-np.linalg.norm(y2))< eps) and (abs(y4n-np.linalg.norm(y4))<eps):
            print('BRERAKING SO THE BOUNDARY CONDITION IS SATISFIED')
            break
        else:
            print('It did not match the boundary condition at infinity!!!')
 
        y2n = np.linalg.norm(y2)
        y4n = np.linalg.norm(y4)
 
    U = y2
    T = y4
    SkinFric = y3
    # Integration for η --> y transformation
    y = np.zeros(N);
    for i in range(N):
       y[i] = y[i-1] + y4[i]*(eta[i]-eta[i-1])
    y = y*np.sqrt(2)
 

    return eta,y, U,T,N, SkinFric 
 
Minf= 0.2
etatest, ytest, Utest, Ttest, Ntest, SkinFric = selfsimilar(Minf)
 
#%%%
plt.plot(Utest[1:],etatest, color='red', label=r'U'.format(Minf))
plt.xlabel("U = df/dEta")
plt.ylabel("Eta")
plt.title('Mach' + ''+ str(Minf))
plt.legend()
plt.show()
plt.grid()

