
#Modelling Nucleus Accumbens
#A Computational Model from Single Cell to Circuit Level

# There are cortex and NAcc groups. The cortex has pyramid and interneurons.
# The NAcc has D1 , D2 type MSNs and interneurons.







from brian2 import *
start_scope()



# Parameters
#neuron parameters
Vthr = 30 * mvolt
EL = -75 * mV




#synaptic parameters
tau_s = 1 * ms
we = 0.1 *amp/mV
wi = 0.1 *amp/mV
Vi = -90 * mV
Ve = 0 * mV
dly=(3+rand())*ms


par_percent=10


#Regular Spike (RS) parameters
#a = 0.02 / ms
#b = 0.2 / ms  
#c = -65 * mvolt        
#d = 8 * mvolt 

#Fast Spike (FS) parameters
#a = 0.1 / ms
#b = 0.2 / ms  
#c = -65 * mvolt        
#d = 2 * mvolt 



#Chattering (CH) parameters
#a = 0.02 / ms
#b = 0.2 / ms  
#c = -50 * mvolt        
#d = 2 * mvolt 




#number_of_neurons_in_cortex=900
number_of_neurons_in_pyramid = 900
number_of_neurons_in_crtx_in = 100

number_of_neurons_in_nacc=450
number_of_neurons_in_msnd1_core = 100
number_of_neurons_in_msnd2_core = 100
number_of_neurons_in_msnd1_shell = 100
number_of_neurons_in_msnd2_shell = 100
number_of_neurons_in_nacc_in = 50

number_of_neurons_in_bg=300
number_of_neurons_in_GPe=100

number_of_neurons_in_thl = 100
number_of_neurons_in_vta = 100



print('Equations')


PG_pyramid = PoissonGroup(number_of_neurons_in_pyramid, 5 * Hz)  #nominal value: 5*Hz

PG_THL = PoissonGroup(number_of_neurons_in_thl, 5 * Hz)

PG_BG = PoissonGroup(number_of_neurons_in_bg, 1 * Hz)


PG_crtx_in = PoissonGroup(number_of_neurons_in_crtx_in, 50 * Hz)    #nominal value: 50*Hz
PG_str_in = PoissonGroup(number_of_neurons_in_nacc_in, 50 * Hz)



eqs_dyn = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v-u)/ms                                           : volt
I : amp
Is=ge*(Ve-v)+gi*(Vi-v) :  amp
dge/dt=-ge/tau_e	: amp/volt
dgi/dt=-gi/tau_i	: amp/volt
a : 1
b : 1
c : volt
d : volt
tau_e : second
tau_i : second
"""


tau_glu=2*ms
tau_DA=1.5*ms
tau_i_msn=tau_s

eqs_msn = """
dv/dt=(0.04/ms/mV)*v**2+(5/ms)*v+140*mV/ms-u/ms+I*mV/(amp*ms)+Is*mV/(amp*ms) : volt
du/dt=a*(b*v+k*mV-u)/ms                                           : volt
I : amp
I_Glu=g_glu*(Ve-v) :  amp
I_DA=g_DA*(V_DA-v) :  amp
I_Ach=g_Ach*(Vi-v)  :  amp
I_GABA=g_GABA*(Vi-v) :  amp
Is=I_Glu+I_DA+I_Ach+I_GABA :  amp
dg_glu/dt=-g_glu/tau_glu	: amp/volt
dg_DA/dt=-g_DA/tau_DA	: amp/volt
dg_Ach/dt=-g_Ach/tau_i_msn	: amp/volt
dg_GABA/dt=-g_GABA/tau_i_msn	: amp/volt
a : 1
b : 1
c : volt
d : volt
k : 1
V_DA : volt
"""




eqs_reset = '''
v = c
u = u+d
'''



pyramid = NeuronGroup(number_of_neurons_in_pyramid, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_pyramid):
    pyramid.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    pyramid.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    pyramid.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    pyramid.d[i] = 8*(100-par_percent+2*par_percent*rand())/100* mvolt 
    pyramid.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    pyramid.u[i] = (-14.5*((100-par_percent+2*par_percent*rand())/100))*mvolt
pyramid.tau_e = tau_s
pyramid.tau_i = tau_s
    
    


crtx_in = NeuronGroup(number_of_neurons_in_crtx_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_crtx_in):
    crtx_in.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    crtx_in.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    crtx_in.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    crtx_in.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    crtx_in.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    crtx_in.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
crtx_in.tau_e = tau_s
crtx_in.tau_i = tau_s

    







msnd1_core =  NeuronGroup(number_of_neurons_in_msnd1_core, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_msnd1_core):
    msnd1_core.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core.c[i] = -62*((100-par_percent+2*par_percent*rand())/100) * mvolt         # NV:-52  
    msnd1_core.d[i] = 0.6*((100-par_percent+2*par_percent*rand())/100) * mvolt          # NV: 1.9
    msnd1_core.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd1_core.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd1_core.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)   #parametre araligi: 0.02 - 0.07
msnd1_core.V_DA = 0*mV



msnd2_core= NeuronGroup(number_of_neurons_in_msnd2_core, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_msnd2_core):
    msnd2_core.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core.c[i] = -60*((100-par_percent+2*par_percent*rand())/100) * mvolt        # NV:-52
    msnd2_core.d[i] = 0.6*((100-par_percent+2*par_percent*rand())/100) * mvolt         #NV: 1.9
    msnd2_core.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd2_core.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd2_core.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
msnd2_core.V_DA = -90*mV




msnd1_shell= NeuronGroup(number_of_neurons_in_msnd1_shell, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd1_shell):
    msnd1_shell.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell.c[i] = -56*((100-par_percent+2*par_percent*rand())/100) * mvolt      # NV:-52  
    msnd1_shell.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt     #NV: 1.9
    msnd1_shell.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd1_shell.u[i] = 35*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd1_shell.k[i] = 35*((100-2*par_percent+4*par_percent*rand())/100)
msnd1_shell.V_DA = 0*mV    
    


msnd2_shell= NeuronGroup(number_of_neurons_in_msnd2_shell, model=eqs_msn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_msnd2_shell):
    msnd2_shell.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell.c[i] = -55*((100-par_percent+2*par_percent*rand())/100) * mvolt        #NV: -52
    msnd2_shell.d[i] = 0.4*((100-par_percent+2*par_percent*rand())/100) * mvolt       #NV: 1.9   
    msnd2_shell.v[i] = (EL-15*mV)*((100-par_percent+2*par_percent*rand())/100)
    msnd2_shell.u[i] = 25*((100-par_percent+2*par_percent*rand())/100)*mV
    msnd2_shell.k[i] = 20*((100-2*par_percent+4*par_percent*rand())/100)
msnd2_shell.V_DA = -90*mV



str_in= NeuronGroup(number_of_neurons_in_nacc_in, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_nacc_in):
    str_in.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    str_in.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    str_in.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    str_in.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt
    str_in.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    str_in.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100)*mV
str_in.tau_e = tau_s
str_in.tau_i = tau_s






GPe = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_GPe):
    GPe.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    GPe.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    GPe.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    GPe.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    GPe.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    GPe.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
GPe.tau_e = tau_s
GPe.tau_i = tau_s*2




GPi = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_GPe):
    GPi.a[i] = 0.01*((100-par_percent+2*par_percent*rand())/100)
    GPi.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    GPi.c[i] = -65*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    GPi.d[i] = 2*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    GPi.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    GPi.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100)*mV
GPi.tau_e = tau_s
GPi.tau_i = tau_s*2





STN = NeuronGroup(number_of_neurons_in_GPe, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)

for i in range(number_of_neurons_in_GPe):
    STN.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    STN.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    STN.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    STN.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    STN.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    STN.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
STN.tau_e = tau_s
STN.tau_i = tau_s







vta_dopamine = NeuronGroup(number_of_neurons_in_vta, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)


for i in range(number_of_neurons_in_vta):
    vta_dopamine.a[i] = 0.02*((100-par_percent+2*par_percent*rand())/100)
    vta_dopamine.b[i] = 0.2*((100-par_percent+2*par_percent*rand())/100)
    vta_dopamine.c[i] = -70*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    vta_dopamine.d[i] = 8*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    vta_dopamine.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    vta_dopamine.u[i] = -14.5 *((100-par_percent+2*par_percent*rand())/100) *mV
vta_dopamine.tau_e = tau_s
vta_dopamine.tau_i = tau_s





##      0.03      0.25    -60     4        0;...      % rebound spike
##      0.03      0.25    -52     0        0;...      % rebound burst
## THL icin rebaund burst alindi.

thl = NeuronGroup(number_of_neurons_in_thl, model=eqs_dyn, method='rk4', threshold='v>Vthr', reset=eqs_reset)



##rebound burst
for i in range(number_of_neurons_in_thl):
    thl.a[i] = 0.03*((100-par_percent+2*par_percent*rand())/100)
    thl.b[i] = 0.25*((100-par_percent+2*par_percent*rand())/100)
    thl.c[i] = -52*((100-par_percent+2*par_percent*rand())/100) * mvolt        
    thl.d[i] = 0.01*((100-par_percent+2*par_percent*rand())/100) * mvolt 
    thl.v[i] = EL*((100-par_percent+2*par_percent*rand())/100)
    thl.u[i] = -14.5*((100-par_percent+2*par_percent*rand())/100) *mV

#regular spike
#for i in range(number_of_neurons_in_thl):
#    thl.a[i]=(0.019*(1+0.1*rand()))
#    thl.b[i]=(0.19*(1+0.1*rand()))
#    thl.c[i]=(-70*(1+0.1*rand())) * mvolt        
#    thl.d[i]=(7.2*(1+0.1*rand())) * mvolt
#    thl.v[i] = EL+1*rand()*mV
#    thl.u[i] = -14.5 *mV+.2*rand()*mV
thl.tau_e = tau_s
thl.tau_i = tau_s*5



        

    


pyramid.I = 0*amp
crtx_in.I = 0*amp

msnd1_core.I = 0*amp
msnd2_core.I = 0*amp
msnd1_shell.I = 0*amp
msnd2_shell.I = 0*amp
str_in.I = 0*amp

GPe.I = 0*amp
GPi.I = 0*amp
STN.I = 0*amp

vta_dopamine.I = 0*amp
thl.I = 0*amp






###### Synapses #############
#############################
print('Synapses')


# 1  cortex
## 11 Pyramid


SP1111 = Synapses(PG_pyramid, pyramid,  'w :siemens', delay=dly, on_pre='ge += w')
SP1111.connect(True, p = 0.25)
SP1111.w=we

S1211 = Synapses(crtx_in, pyramid,  delay=dly, on_pre='gi += wi')
S1211.connect(True, p = 0.25)



# 12 IN

SP1212 = Synapses(PG_crtx_in, crtx_in,  delay=dly, on_pre='ge += 5*we')
SP1212.connect(True, p = 0.25)

S1112 = Synapses(pyramid, crtx_in, delay=dly, on_pre='ge += we')
S1112.connect(True, p = 0.25)





# 2 NAcc
# 21 Core MSND1

S1121 = Synapses(pyramid, msnd1_core, 'w :siemens', delay=dly, on_pre='g_glu += w')
S1121.connect(True, p = 0.25)
S1121.w=we
print("synapse number from cortex: "+str(S1121.N)+" average:"+str(S1121.N/100))

S2521 = Synapses(str_in, msnd1_core, 'w :siemens', delay=dly, on_pre='g_Ach += w')
S2521.connect(True, p = 0.25)
S2521.w=wi
print("synapse number from IN: "+str(S2521.N)+" average:"+str(S2521.N/100))

S4121 = Synapses( thl, msnd1_core, 'w :siemens', delay=dly, on_pre='g_glu += w')
S4121.connect(True, p = 0.25)
S4121.w=we
print("synapse number from THL: "+str(S4121.N)+" average:"+str(S4121.N/100))

S5121 = Synapses(vta_dopamine, msnd1_core, 'w :siemens', delay=dly, on_pre='g_DA += w')
S5121.connect(True, p = 0.25)
S5121.w=we
print("synapse number from VTA: "+str(S5121.N)+" average:"+str(S5121.N/100))

##### Colateral inhibitions from MSNs

S2121 = Synapses(msnd1_core, msnd1_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2121.connect(True, p = 0.05)
S2121.w=wi
print("synapse number from MSND1C: "+str(S2121.N)+" average:"+str(S2121.N/100))

S2221 = Synapses(msnd2_core, msnd1_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2221.connect(True, p = 0.25)
S2221.w=wi*2
print("synapse number from MSND2C: "+str(S2221.N)+" average:"+str(S2221.N/100))

S2321 = Synapses(msnd1_shell, msnd1_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2321.connect(True, p = 0.05)
S2321.w=wi
print("synapse number from MSND1S: "+str(S2321.N)+" average:"+str(S2321.N/100))

S2421 = Synapses(msnd2_shell, msnd1_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2421.connect(True, p = 0.25)
S2421.w=wi*2
print("synapse number from MSND2S: "+str(S2421.N)+" average:"+str(S2421.N/100))
print("Total average synapse number:" +str((S1121.N+S2521.N+S4121.N+S5121.N+S2121.N+S2221.N+S2321.N+S2421.N)/100))





# 22 Core MSND2

S1122 = Synapses(pyramid, msnd2_core,  'w :siemens', delay=dly, on_pre='g_glu += w')
S1122.connect(True, p = 0.25)
S1122.w=we

S2522 = Synapses(str_in, msnd2_core,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S2522.connect(True, p = 0.25)
S2522.w=wi

S4122 = Synapses( thl, msnd2_core,  'w :siemens', delay=dly, on_pre='g_glu += w')
S4122.connect(True, p = 0.25)
S4122.w=we

S5122 = Synapses(vta_dopamine, msnd2_core,  'w :siemens', delay=dly, on_pre='g_DA += w')
S5122.connect(True, p = 0.25)
S5122.w=wi


##### Colateral inhibitions from MSNs

S2122 = Synapses(msnd1_core, msnd2_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2122.connect(True, p = 0.25)
S2122.w=wi*2


S2222 = Synapses(msnd2_core, msnd2_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2222.connect(True, p = 0.05)
S2222.w=wi


S2322 = Synapses(msnd1_shell, msnd2_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2322.connect(True, p = 0.25)
S2322.w=wi*2


S2422 = Synapses(msnd2_shell, msnd2_core, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2422.connect(True, p = 0.05)
S2422.w=wi





# 23 Shell MSND1

S1123 = Synapses(pyramid, msnd1_shell,  'w :siemens', delay=dly, on_pre='g_glu += w')
S1123.connect(True, p = 0.25)
S1123.w=we

S2523 = Synapses(str_in, msnd1_shell,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S2523.connect(True, p = 0.25)
S2523.w=wi

S4123 = Synapses(thl, msnd1_shell,  'w :siemens', delay=dly, on_pre='g_glu += w')
S4123.connect(True, p = 0.25)
S4123.w=we

S5123 = Synapses(vta_dopamine, msnd1_shell,  'w :siemens', delay=dly, on_pre='g_DA += w')
S5123.connect(True, p = 0.25)
S5123.w=we




##### Colateral inhibitions from MSNs

S2123 = Synapses(msnd1_core, msnd1_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2123.connect(True, p = 0.05)
S2123.w=wi


S2223 = Synapses(msnd2_core, msnd1_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2223.connect(True, p = 0.25)
S2223.w=wi*2


S2323 = Synapses(msnd1_shell, msnd1_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2323.connect(True, p = 0.05)
S2323.w=wi


S2423 = Synapses(msnd2_shell, msnd1_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2423.connect(True, p = 0.25)
S2423.w=wi*2






# 24 Shell MSND2

S1124 = Synapses(pyramid, msnd2_shell,  'w :siemens', delay=dly, on_pre='g_glu += w')
S1124.connect(True, p = 0.25)
S1124.w=we

S2524 = Synapses(str_in, msnd2_shell,  'w :siemens', delay=dly, on_pre='g_Ach += w')
S2524.connect(True, p = 0.25)
S2524.w=wi

S4124 = Synapses(thl, msnd2_shell,  'w :siemens', delay=dly, on_pre='g_glu += w')
S4124.connect(True, p = 0.25)
S4124.w=we

S5124 = Synapses(vta_dopamine, msnd2_shell,  'w :siemens', delay=dly, on_pre='g_DA += w')
S5124.connect(True, p = 0.25)
S5124.w=wi



##### Colateral inhibitions from MSNs

S2124 = Synapses(msnd1_core, msnd2_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2124.connect(True, p = 0.25)
S2124.w=wi*2


S2224 = Synapses(msnd2_core, msnd2_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2224.connect(True, p = 0.05)
S2224.w=wi


S2324 = Synapses(msnd1_shell, msnd2_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2324.connect(True, p = 0.25)
S2324.w=wi*2


S2424 = Synapses(msnd2_shell, msnd2_shell, 'w :siemens', delay=dly, on_pre='g_GABA += w')
S2424.connect(True, p = 0.05)
S2424.w=wi





# 25 IN


S1125 = Synapses(pyramid, str_in,  delay=dly, on_pre='ge += 0.25*we')
S1125.connect(True, p = 0.20)

S2125 = Synapses(msnd1_core, str_in, delay=dly, on_pre='gi += wi')
S2125.connect(True, p = 0.25)

S2225 = Synapses(msnd2_core, str_in, delay=dly, on_pre='gi += wi')
S2225.connect(True, p = 0.25)

S2325 = Synapses(msnd1_shell, str_in, delay=dly, on_pre='gi += wi')
S2325.connect(True, p = 0.25)

S2425 = Synapses(msnd2_shell, str_in, delay=dly, on_pre='gi += wi')
S2425.connect(True, p = 0.25)





 




# 31 GPe

SP3031 = Synapses(PG_BG, GPe,  delay=dly, on_pre='ge += 2*we')
SP3031.connect(True, p = 0.25)

S2231 = Synapses(msnd2_core, GPe,  delay=dly, on_pre='gi += wi')
S2231.connect(True, p = 0.25)

S2431 = Synapses(msnd2_shell, GPe,  delay=dly, on_pre='gi += wi')
S2431.connect(True, p = 0.25)

S3331 = Synapses(STN, GPe,  delay=dly, on_pre='ge += we')
S3331.connect(True, p = 0.25)




# 32 GPi

SP3032 = Synapses(PG_BG, GPi,  delay=dly, on_pre='ge += 2*we')
SP3032.connect(True, p = 0.25)

S2132 = Synapses(msnd1_core, GPi,  delay=dly, on_pre='gi += wi')
S2132.connect(True, p = 0.25)

S2332 = Synapses(msnd1_shell, GPi,  delay=dly, on_pre='gi += wi')
S2332.connect(True, p = 0.25)

S3132 = Synapses(GPe, GPi,  delay=dly, on_pre='gi += wi')
S3132.connect(True, p = 0.25)







# 33 STN

S1133 = Synapses(pyramid, STN,  delay=dly, on_pre='ge += 0.5*we')
S1133.connect(True, p = 0.25)

S3133 = Synapses(GPe, STN,  delay=dly, on_pre='gi += wi')
S3133.connect(True, p = 0.25)





# 41 THL

SP4141 = Synapses(PG_THL, thl,  delay=dly, on_pre='ge += 0.5*we')
SP4141.connect(True, p = 0.25)

S3241 = Synapses(GPi, thl, delay=dly, on_pre='gi += 0.5*wi')
S3241.connect(True, p = 0.25)






# 51 VTA

SP5151 = Synapses(PG_pyramid, vta_dopamine,  'w :siemens', delay=dly, on_pre='ge += w')
SP5151.connect(True, p = 0.25)
SP5151.w=we     



#==============================================================================
#==============================================================================
#==============================================================================






import time
init_time=time.time()


######### First 100ms #########
################################
duration1=100*ms


SP1111.w=we*0
SP5151.w=we*0

print('100 ms initial conditions delay')
print("sim_time="+str(duration1))
run(duration1, report='text')

SP1111.w=we
SP5151.w=we

###### Monitors ###########
###########################
print('Monitors')




trace_pyramid = StateMonitor(pyramid, 'v', record=9)
#monge_pyramid = StateMonitor(pyramid, 'ge', record=True)
#mongi_pyramid = StateMonitor(pyramid, 'gi', record=True)
spikes_pyramid = SpikeMonitor(pyramid)


trace_crtx_in = StateMonitor(crtx_in, 'v', record=9)
#monge_crtx_in = StateMonitor(crtx_in, 'ge', record=True)
#mongi_crtx_in = StateMonitor(crtx_in, 'gi', record=True)
spikes_crtx_in = SpikeMonitor(crtx_in)







trace_msnd1_core = StateMonitor(msnd1_core, 'v', record=True)
trace_g_glu_msnd1_core = StateMonitor(msnd1_core, 'g_glu', record=9)
trace_g_DA_msnd1_core = StateMonitor(msnd1_core, 'g_DA', record=9)
trace_g_Ach_msnd1_core = StateMonitor(msnd1_core, 'g_Ach', record=9)
trace_g_GABA_msnd1_core = StateMonitor(msnd1_core, 'g_GABA', record=9)
trace_I_Glu_msnd1_core = StateMonitor(msnd1_core, 'I_Glu', record=9)
trace_I_DA_msnd1_core = StateMonitor(msnd1_core, 'I_DA', record=9)
trace_I_Ach_msnd1_core = StateMonitor(msnd1_core, 'I_Ach', record=9)
trace_I_GABA_msnd1_core = StateMonitor(msnd1_core, 'I_GABA', record=9)
trace_I_s_msnd1_core = StateMonitor(msnd1_core, 'Is', record=True)
spikes_msnd1_core = SpikeMonitor(msnd1_core)






trace_msnd2_core = StateMonitor(msnd2_core, 'v', record=True)
trace_g_glu_msnd2_core = StateMonitor(msnd2_core, 'g_glu', record=9)
trace_g_DA_msnd2_core = StateMonitor(msnd2_core, 'g_DA', record=9)
trace_g_Ach_msnd2_core = StateMonitor(msnd2_core, 'g_Ach', record=9)
trace_g_GABA_msnd2_core = StateMonitor(msnd2_core, 'g_GABA', record=9)
trace_I_Glu_msnd2_core = StateMonitor(msnd2_core, 'I_Glu', record=9)
trace_I_DA_msnd2_core = StateMonitor(msnd2_core, 'I_DA', record=9)
trace_I_Ach_msnd2_core = StateMonitor(msnd2_core, 'I_Ach', record=9)
trace_I_GABA_msnd2_core = StateMonitor(msnd2_core, 'I_GABA', record=9)
trace_I_s_msnd2_core = StateMonitor(msnd2_core, 'Is', record=True)
spikes_msnd2_core = SpikeMonitor(msnd2_core)


trace_msnd1_shell = StateMonitor(msnd1_shell, 'v', record=True)
#monge_msnd1_shell = StateMonitor(msnd1_shell, 'ge', record=True)
#mongi_msnd1_shell = StateMonitor(msnd1_shell, 'gi', record=True)
spikes_msnd1_shell = SpikeMonitor(msnd1_shell)
trace_I_s_msnd1_shell = StateMonitor(msnd1_shell, 'Is', record=True)


trace_msnd2_shell = StateMonitor(msnd2_shell, 'v', record=True)
#monge_msnd2_shell = StateMonitor(msnd2_shell, 'ge', record=True)
#mongi_msnd2_shell = StateMonitor(msnd2_shell, 'gi', record=True)
spikes_msnd2_shell = SpikeMonitor(msnd2_shell)
trace_I_s_msnd2_shell = StateMonitor(msnd2_shell, 'Is', record=True)


trace_str_in = StateMonitor(str_in, 'v', record=True)
#monge_str_in = StateMonitor(str_in, 'ge', record=True)
#mongi_str_in = StateMonitor(str_in, 'gi', record=True)
spikes_str_in = SpikeMonitor(str_in)
trace_I_s_nacc_in = StateMonitor(str_in, 'Is', record=True)

trace_GPe = StateMonitor(GPe, 'v', record=9)
#monge_GPe = StateMonitor(GPe, 'ge', record=True)
#mongi_GPe = StateMonitor(GPe, 'gi', record=True)
spikes_GPe = SpikeMonitor(GPe)


trace_GPi = StateMonitor(GPi, 'v', record=9)
#monge_GPi = StateMonitor(GPi, 'ge', record=True)
#mongi_GPi = StateMonitor(GPi, 'gi', record=True)
spikes_GPi = SpikeMonitor(GPi)



trace_STN = StateMonitor(STN, 'v', record=9)
#monge_STN = StateMonitor(STN, 'ge', record=True)
#mongi_STN = StateMonitor(STN, 'gi', record=True)
spikes_STN = SpikeMonitor(STN)




trace_vta_dopamine = StateMonitor(vta_dopamine, 'v', record=9)
#monge_vta_dopamine = StateMonitor(vta_dopamine, 'ge', record=True)
#mongi_vta_dopamine = StateMonitor(vta_dopamine, 'gi', record=True)
spikes_vta_dopamine = SpikeMonitor(vta_dopamine)


trace_thl = StateMonitor(thl, 'v', record=9)
#monge_pyramid = StateMonitor(pyramid, 'ge', record=True)
#mongi_pyramid = StateMonitor(pyramid, 'gi', record=True)
spikes_thl = SpikeMonitor(thl)



########### Excitation currents ##########
##########################################

print("start")


SP1111.w=we*2
SP5151.w=we


number_of_scenario=1




#########################  ----- Scenarios   --------------  ####################

##-------------------- Scenario 0 --------------------------------##
## It is used for testing purposes.

if number_of_scenario==0:

    wcne = 0.20 *we   
    wtne = 0.25 *we   
    wvne = 3.00 *we
    wvni = 3.00 *wi
    
    SP1111.w=we
    SP5151.w=we
        
    
    S1121.w=wcne        # from pyramid to msnd1 core
    S4121.w=wtne        # from THL to     msnd1 core
    S5121.w=wvne        # from VTA to     msnd1 core
    
    
    S1122.w=wcne        # from pyramid to msnd2 core
    S4122.w=wtne        # from THL to     msnd2 core
    S5122.w=wvni        # from VTA to     msnd2 core
    
    
    S1123.w=wcne        # from pyramid to msnd1 shell
    S4123.w=wtne        # from THL to     msnd1 shell
    S5123.w=wvne        # from VTA to     msnd1 shell
    
    
    S1124.w=wcne        # from pyramid to msnd2 shell
    S4124.w=wtne        # from THL to     msnd2 shell
    S5124.w=wvni        # from VTA to     msnd2 shell
    
    
    
    print('Scenario 0')
    run(100*ms,report='text')



##-------------------- Scenario 1 --------------------------------##

##### Changing cortex and VTA input, investigate the relation of stimulus and reward

elif number_of_scenario==1:


    print('Scenario 1')
    duration1=100*ms
    duration2=1000*ms
    
    wcne = 0.20 *we   
    wtne = 0.25 *we   
    wvne = 3.00 *we
    wvni = 3.00 *wi
    
    SP1111.w=we*0
    SP5151.w=we*0
            
    
    S1121.w=wcne        # from pyramid to msnd1 core
    S4121.w=wtne        # from THL to     msnd1 core
    S5121.w=wvne        # from VTA to     msnd1 core
    
    
    S1122.w=wcne        # from pyramid to msnd2 core
    S4122.w=wtne        # from THL to     msnd2 core
    S5122.w=wvni        # from VTA to     msnd2 core
    
    
    S1123.w=wcne        # from pyramid to msnd1 shell
    S4123.w=wtne        # from THL to     msnd1 shell
    S5123.w=wvne        # from VTA to     msnd1 shell
    
    
    S1124.w=wcne        # from pyramid to msnd2 shell
    S4124.w=wtne        # from THL to     msnd2 shell
    S5124.w=wvni        # from VTA to     msnd2 shell
    
    
    print('Scenario 1: D1 Resting')    
    
    S1121.w=wcne*0        # from pyramid to msnd1 core
    S1122.w=wcne*0        # from pyramid to msnd2 core    
    S1123.w=wcne*0        # from pyramid to msnd1 shell    
    S1124.w=wcne*0        # from pyramid to msnd2 shell
    SP1111.w=we*0
    SP5151.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    print('Scenario 1: D2  Only Cortex')    

    S1121.w=wcne        # from pyramid to msnd1 core
    S1122.w=wcne        # from pyramid to msnd2 core    
    S1123.w=wcne        # from pyramid to msnd1 shell    
    S1124.w=wcne        # from pyramid to msnd2 shell    
    SP1111.w=we
    SP5151.w=we*0.1
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')    
    
    
    
    print('Scenario 1: D3 Resting')    
    
    SP1111.w=we*0
    SP5151.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')
    
    
    print('Scenario 1: D4 Only VTA')    
    
    S1121.w=wcne*0.1        # from pyramid to msnd1 core
    S1122.w=wcne*0.1        # from pyramid to msnd2 core    
    S1123.w=wcne*0.1        # from pyramid to msnd1 shell    
    S1124.w=wcne*0.1        # from pyramid to msnd2 shell    
    SP1111.w=we
    SP5151.w=we
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')

    
    print('Scenario 1: D5 Resting')    
    
    SP1111.w=we*0
    SP5151.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')


    print('Scenario 1: D6 Cortex and VTA')    
    
    S1121.w=wcne        # from pyramid to msnd1 core
    S1122.w=wcne        # from pyramid to msnd2 core    
    S1123.w=wcne        # from pyramid to msnd1 shell    
    S1124.w=wcne        # from pyramid to msnd2 shell    
    SP1111.w=we
    SP5151.w=we
    
    print("sim_time="+str(duration2))
    run(duration2, report='text')

    
    print('Scenario 1: D7 Resting')    
    
    SP1111.w=we*0
    SP5151.w=we*0
    
    print("sim_time="+str(duration1))
    run(duration1, report='text')



#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####

##-------------------- Scenario 3 --------------------------------##

 
else:
    print('No Scenario!!!!')


    

final_time=time.time()
print("Simulation Time:",str(final_time-init_time))

#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####
#### ------------------------------------------------------------------- ####


print("Pyramid : "+str(spikes_pyramid.num_spikes)+"  ---->  : " +str(spikes_pyramid.num_spikes/number_of_neurons_in_pyramid))
print("Crtx IN : "+str(spikes_crtx_in.num_spikes)+"  ---->  : " +str(spikes_crtx_in.num_spikes/number_of_neurons_in_crtx_in))
print("MSND1 Core : "+str(spikes_msnd1_core.num_spikes)+"  ---->  : " +str(spikes_msnd1_core.num_spikes/number_of_neurons_in_msnd1_core))
print("MSND2 Core : "+str(spikes_msnd2_core.num_spikes)+"  ---->  : " +str(spikes_msnd2_core.num_spikes/number_of_neurons_in_msnd2_core))
print("MSND1 Shell : "+str(spikes_msnd1_shell.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell.num_spikes/number_of_neurons_in_msnd1_shell))
print("MSND2 Shell : "+str(spikes_msnd2_shell.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell.num_spikes/number_of_neurons_in_msnd2_shell))
print("NAcc IN : "+str(spikes_str_in.num_spikes)+"  ---->  : " +str(spikes_pyramid.num_spikes/number_of_neurons_in_nacc_in))
print("GPe : "+str(spikes_GPe.num_spikes)+"  ---->  : " +str(spikes_GPe.num_spikes/number_of_neurons_in_GPe))
print("GPi : "+str(spikes_GPi.num_spikes)+"  ---->  : " +str(spikes_GPi.num_spikes/number_of_neurons_in_GPe))
print("STN : "+str(spikes_STN.num_spikes)+"  ---->  : " +str(spikes_STN.num_spikes/number_of_neurons_in_GPe))
print("VTA : "+str(spikes_vta_dopamine.num_spikes)+"  ---->  : " +str(spikes_vta_dopamine.num_spikes/number_of_neurons_in_vta))
print("THL : "+str(spikes_thl.num_spikes)+"  ---->  : " +str(spikes_thl.num_spikes/number_of_neurons_in_thl))




#########################################
############# --- Figures --- #########
#########################################

figure()
subplot(411)
plot(trace_pyramid.t / ms, trace_pyramid[9].v / mV)
ylabel('Pyramid')


subplot(412)
plot(trace_crtx_in.t / ms, trace_crtx_in[9].v / mV)
ylabel('Crtx IN')


subplot(413)
plot(trace_vta_dopamine.t / ms, trace_vta_dopamine[9].v / mV)
ylabel('VTA')


subplot(414)
plot(trace_thl.t / ms, trace_thl[9].v / mV, linewidth=1)
xlabel('time (ms)')
ylabel('THL')
tight_layout()
savefig('example1.pdf', dpi=600)

figure()
subplot(411)
plot(trace_msnd1_core.t / ms, trace_msnd1_core[9].v / mV)
ylabel('MSND1 Core')



subplot(412)
plot(trace_msnd2_core.t / ms, trace_msnd2_core[9].v / mV)
ylabel('MSND2 Core')




subplot(413)
plot(trace_msnd1_shell.t / ms, trace_msnd1_shell[9].v / mV)
ylabel('MSND1 Shell')



subplot(414)
plot(trace_msnd2_shell.t / ms, trace_msnd2_shell[9].v / mV)
xlabel('time, ms')
ylabel('MSND2 Shell')






figure()



subplot(511)
plot(trace_I_Glu_msnd1_core.t / ms,trace_I_Glu_msnd1_core[9].I_Glu, label='monitorden')
ylabel('gGlu')



subplot(512)
plot(trace_I_DA_msnd1_core.t / ms,trace_I_DA_msnd1_core[9].I_DA, label='monitorden')
ylabel('gDA')




subplot(513)
plot(trace_I_GABA_msnd1_core.t / ms,trace_I_GABA_msnd1_core[9].I_GABA, label='monitorden')
ylabel('gGABA')



subplot(514)
plot(trace_I_Ach_msnd1_core.t / ms,trace_I_Ach_msnd1_core[9].I_Ach, label='monitorden')
ylabel('gAch')



subplot(515)
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[9].Is)
hlines(0,100,300)
xlabel('time, ms')
ylabel('Is')


figure()
plot(trace_msnd1_core.t / ms, trace_msnd1_core[9].v / mV, label='v')
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[9].Is, label='Is')
plot(trace_I_Glu_msnd1_core.t / ms,trace_I_Glu_msnd1_core[9].I_Glu, label='Glu')
plot(trace_I_DA_msnd1_core.t / ms,trace_I_DA_msnd1_core[9].I_DA,label='DA')
plot(trace_I_GABA_msnd1_core.t / ms,trace_I_GABA_msnd1_core[9].I_GABA, label='GABA')
plot(trace_I_Ach_msnd1_core.t / ms,trace_I_Ach_msnd1_core[9].I_Ach,label='Ach')




figure()

subplot(511)
plot(trace_I_Glu_msnd2_core.t / ms,trace_I_Glu_msnd2_core[9].I_Glu, label='monitorden')
ylabel('gGlu2')



subplot(512)
plot(trace_I_DA_msnd2_core.t / ms,trace_I_DA_msnd2_core[9].I_DA, label='monitorden')
ylabel('gDA2')




subplot(513)
plot(trace_I_GABA_msnd2_core.t / ms,trace_I_GABA_msnd2_core[9].I_GABA, label='monitorden')
ylabel('gGABA2')



subplot(514)
plot(trace_I_Ach_msnd2_core.t / ms,trace_I_Ach_msnd2_core[9].I_Ach, label='monitorden')
ylabel('gAch2')



subplot(515)
plot(trace_I_s_msnd2_core.t / ms, trace_I_s_msnd2_core[9].Is)
hlines(0,100,300)
xlabel('time, ms')
ylabel('Is2')


figure()
plot(trace_msnd2_core.t / ms, trace_msnd2_core[9].v / mV,label='v')
plot(trace_I_s_msnd2_core.t / ms, trace_I_s_msnd2_core[9].Is,label='Is')
plot(trace_I_Glu_msnd2_core.t / ms,trace_I_Glu_msnd2_core[9].I_Glu,label='Glu')
plot(trace_I_DA_msnd2_core.t / ms,trace_I_DA_msnd2_core[9].I_DA, label='DA')
plot(trace_I_GABA_msnd2_core.t / ms,trace_I_GABA_msnd2_core[9].I_GABA, label='GABA')
plot(trace_I_Ach_msnd2_core.t / ms,trace_I_Ach_msnd2_core[9].I_Ach, label='Ach')



figure()
subplot(411)
plot(trace_GPe.t / ms, trace_GPe[9].v / mV)
ylabel('GPe')

subplot(412)
plot(trace_GPi.t / ms, trace_GPi[9].v / mV)
ylabel('GPi')

subplot(413)
plot(trace_STN.t / ms, trace_STN[9].v / mV)
ylabel('STN')

subplot(414)
plot(trace_str_in.t / ms, trace_str_in[9].v / mV)
xlabel('time (ms)')
ylabel('NAcc IN')


#==============================================================================


figure()
subplot(411)
plot(spikes_pyramid.t/ms, spikes_pyramid.i, '.k')

ylabel('Pyramid');

subplot(412)
plot(spikes_crtx_in.t/ms, spikes_crtx_in.i, '.k')

ylabel('Crtx IN');


subplot(413)
plot(spikes_vta_dopamine.t/ms, spikes_vta_dopamine.i, '.k')
ylabel('VTA');


subplot(414)
plot(spikes_thl.t/ms, spikes_thl.i, '.k')
xlabel('Time (ms)')
ylabel('THL');



figure()

subplot(511)
plot(spikes_msnd1_core.t/ms, spikes_msnd1_core.i, '.k')
ylabel('MSND1 Core');

subplot(512)
plot(spikes_msnd2_core.t/ms, spikes_msnd2_core.i, '.k')
ylabel('MSND2 Core');


subplot(513)
plot(spikes_msnd1_shell.t/ms, spikes_msnd1_shell.i, '.k')
ylabel('MSND1 Shell');

subplot(514)
plot(spikes_msnd2_shell.t/ms, spikes_msnd2_shell.i, '.k')
ylabel('MSND2 Shell');


subplot(515)
plot(spikes_str_in.t/ms, spikes_str_in.i, '.k')
xlabel('time, ms')
ylabel('NAcc IN');



figure()
subplot(311)
plot(spikes_GPe.t/ms, spikes_GPe.i, '.k')
ylabel('GPe');

subplot(312)
plot(spikes_GPi.t/ms, spikes_GPi.i, '.k')
ylabel('GPi');


subplot(313)
plot(spikes_STN.t/ms, spikes_STN.i, '.k')
ylabel('STN');


################################################################
################################################################
#

figure()
subplot(421)
plot(spikes_pyramid.t/ms, spikes_pyramid.i, '.k')

ylabel('Pyramid');



subplot(422)
#plot(spikes_crtx_in.t/ms, spikes_crtx_in.i, '.k')
plot(spikes_vta_dopamine.t/ms, spikes_vta_dopamine.i, '.k')
ylabel('VTA');


subplot(423)
plot(spikes_msnd1_core.t/ms, spikes_msnd1_core.i, '.k')
ylabel('MSND1 Core');

subplot(424)
plot(spikes_msnd2_core.t/ms, spikes_msnd2_core.i, '.k')
ylabel('MSND2 Core');


subplot(425)
plot(spikes_msnd1_shell.t/ms, spikes_msnd1_shell.i, '.k')
ylabel('MSND1 Shell');

subplot(426)
plot(spikes_msnd2_shell.t/ms, spikes_msnd2_shell.i, '.k')
ylabel('MSND2 Shell');

subplot(427)
plot(spikes_STN.t/ms, spikes_STN.i, '.k')
ylabel('STN');

subplot(428)
plot(spikes_thl.t/ms, spikes_thl.i, '.k')
xlabel('Time (ms)')
ylabel('THL');


figure()

subplot(511)
hist_cortex = hist(spikes_pyramid.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
ylabel('Cortex fr (sp/s)');

subplot(512)
hist_vta_dopamine = hist(spikes_vta_dopamine.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
ylabel('vta fr (sp/s)');

subplot(513)

hist_msnd1_core = hist(spikes_msnd1_core.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
hist_msnd1_shell = hist(spikes_msnd1_shell.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
ylabel('msnd1');
#
subplot(514)
hist_msnd2_core = hist(spikes_msnd2_core.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
hist_msnd2_shell = hist(spikes_msnd2_shell.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
ylabel('msnd2');


subplot(515)
hist_thl = hist(spikes_thl.t/ms, 100)#, histtype='stepfilled', facecolor='k', weights=ones(len(spikes_str_in))/(number_of_neurons_in_str_in*defaultclock.dt))
ylabel('thl fr (sp/s)');

show()


#==============================================================================
#========     Firing Rate data ====================================================
#==============================================================================


fr_cortex=hist_cortex[0]


fr_msnd1c=hist_msnd1_core[0]
fr_msnd2c=hist_msnd2_core[0]
fr_msnd1s=hist_msnd1_shell[0]
fr_msnd2s=hist_msnd2_shell[0]


fr_vta_dopamine=hist_vta_dopamine[0]

fr_thl=hist_thl[0]
#==============================================================================
# fr_nacc_in=hist_str_in[0]
#==============================================================================


fr_nacc=fr_msnd1c+fr_msnd2c+fr_msnd1s+fr_msnd2s
fr_D1=fr_msnd1c+fr_msnd1s
fr_D2=fr_msnd2c+fr_msnd2s

#==============================================================================
#========     File Writing ====================================================
#==============================================================================

filename=str(time.time())


filename_fundamentals_informations="001_filename_fundamentals_informations_"+filename+".dat"
with open(filename_fundamentals_informations,"w") as results_filename_fundamentals_informations:
    results_filename_fundamentals_informations.write("benzetim suresi : "+str(final_time-init_time)+"\n")
    results_filename_fundamentals_informations.write("Pyramid : "+str(spikes_pyramid.num_spikes)+"  ---->  : " +str(spikes_pyramid.num_spikes/number_of_neurons_in_pyramid)+"\n")
    results_filename_fundamentals_informations.write("Crtx IN : "+str(spikes_crtx_in.num_spikes)+"  ---->  : " +str(spikes_crtx_in.num_spikes/number_of_neurons_in_crtx_in)+"\n")
    results_filename_fundamentals_informations.write("MSND1 Core : "+str(spikes_msnd1_core.num_spikes)+"  ---->  : " +str(spikes_msnd1_core.num_spikes/number_of_neurons_in_msnd1_core)+"\n")
    results_filename_fundamentals_informations.write("MSND2 Core : "+str(spikes_msnd2_core.num_spikes)+"  ---->  : " +str(spikes_msnd2_core.num_spikes/number_of_neurons_in_msnd2_core)+"\n")
    results_filename_fundamentals_informations.write("MSND1 Shell : "+str(spikes_msnd1_shell.num_spikes)+"  ---->  : " +str(spikes_msnd1_shell.num_spikes/number_of_neurons_in_msnd1_shell)+"\n")
    results_filename_fundamentals_informations.write("MSND2 Shell : "+str(spikes_msnd2_shell.num_spikes)+"  ---->  : " +str(spikes_msnd2_shell.num_spikes/number_of_neurons_in_msnd2_shell)+"\n")
    results_filename_fundamentals_informations.write("NAcc IN : "+str(spikes_str_in.num_spikes)+"  ---->  : " +str(spikes_pyramid.num_spikes/number_of_neurons_in_nacc_in)+"\n")
    results_filename_fundamentals_informations.write("GPe : "+str(spikes_GPe.num_spikes)+"  ---->  : " +str(spikes_GPe.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("GPi : "+str(spikes_GPi.num_spikes)+"  ---->  : " +str(spikes_GPi.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("STN : "+str(spikes_STN.num_spikes)+"  ---->  : " +str(spikes_STN.num_spikes/number_of_neurons_in_GPe)+"\n")
    results_filename_fundamentals_informations.write("VTA : "+str(spikes_vta_dopamine.num_spikes)+"  ---->  : " +str(spikes_vta_dopamine.num_spikes/number_of_neurons_in_vta)+"\n")
    results_filename_fundamentals_informations.write("THL : "+str(spikes_thl.num_spikes)+"  ---->  : " +str(spikes_thl.num_spikes/number_of_neurons_in_thl)+"\n")




filename_crtx="011_crtx_"+filename+".dat"
with open(filename_crtx,"w") as results_crtx:
    for i in range(len(fr_cortex)):
 		 	results_crtx.write(str(fr_cortex[i])+"\n")


filename_str="012_str_"+filename+".dat"
with open(filename_str,"w") as results_str:
    for i in range(len(fr_nacc)):
 		 	results_str.write(str(fr_nacc[i])+"\n")
     
     

filename_msnd1c="001_trace_msnd1_core_"+filename+".dat"
with open(filename_msnd1c,"w") as results_msnd1c:
    for i in range(len(trace_msnd1_core.t)):
 		 	results_msnd1c.write(str(trace_msnd1_core.t[i]/ms)+","+str(trace_msnd1_core[9].v[i]/mV)+"\n")


filename_msnd2c="002_trace_msnd2_core_"+filename+".dat"
with open(filename_msnd2c,"w") as results_msnd2c:
    for i in range(len(trace_msnd2_core.t)):
 		 	results_msnd2c.write(str(trace_msnd2_core.t[i]/ms)+","+str(trace_msnd2_core[9].v[i]/mV)+"\n")

     
filename_msnd1s="003_trace_msnd1_shell_"+filename+".dat"
with open(filename_msnd1s,"w") as results_msnd1s:
    for i in range(len(trace_msnd1_shell.t)):
 		 	results_msnd1s.write(str(trace_msnd1_shell.t[i]/ms)+","+str(trace_msnd1_shell[9].v[i]/mV)+"\n")
              
filename_msnd2s="004_trace_msnd2_shell_"+filename+".dat"
with open(filename_msnd2s,"w") as results_msnd2s:
    for i in range(len(trace_msnd2_shell.t)):
 		 	results_msnd2s.write(str(trace_msnd2_shell.t[i]/ms)+","+str(trace_msnd2_shell[9].v[i]/mV)+"\n")              
     
     

filename_innacc="005_trace_in_nacc_"+filename+".dat"
with open(filename_innacc,"w") as results_innacc:
    for i in range(len(trace_str_in.t)):
 		 	results_innacc.write(str(trace_str_in.t[i]/ms)+","+str(trace_str_in[9].v[i]/mV)+"\n")


filename_msnd1c_Is="001_trace_I_s_msnd1_core"+filename+".dat"
with open(filename_msnd1c_Is,"w") as results_msnd1c_Is:
    for i in range(len(trace_I_s_msnd1_core.t)):
 		 	results_msnd1c_Is.write(str(trace_I_s_msnd1_core.t[i]/ms)+","+str(trace_I_s_msnd1_core[9].Is[i]/amp)+","+str(trace_I_Glu_msnd1_core[9].I_Glu[i]/amp)+","+str(trace_I_DA_msnd1_core[9].I_DA[i]/amp)+","+str(trace_I_GABA_msnd1_core[9].I_GABA[i]/amp)+","+str(trace_I_Ach_msnd1_core[9].I_Ach[i]/amp)+"\n")


filename_msnd2c_Is="002_trace_I_s_msnd2_core"+filename+".dat"
with open(filename_msnd2c_Is,"w") as results_msnd2c_Is:
    for i in range(len(trace_I_s_msnd2_core.t)):
 		 	results_msnd2c_Is.write(str(trace_I_s_msnd2_core.t[i]/ms)+","+str(trace_I_s_msnd2_core[9].Is[i]/amp)+","+str(trace_I_Glu_msnd2_core[9].I_Glu[i]/amp)+","+str(trace_I_DA_msnd2_core[9].I_DA[i]/amp)+","+str(trace_I_GABA_msnd2_core[9].I_GABA[i]/amp)+","+str(trace_I_Ach_msnd2_core[9].I_Ach[i]/amp)+"\n")

filename_rasterplot_msnd1c="001_rasterplot_msnd1_core_"+filename+".dat"
with open(filename_rasterplot_msnd1c,"w") as results_rasterplot__msnd1c:
    for i_c in range(len(spikes_msnd1_core.t/ms)):
        results_rasterplot__msnd1c.write(str(spikes_msnd1_core.t[i_c]/ms)+","+str(spikes_msnd1_core.i[i_c])+"\n")

        
filename_rasterplot_msnd2c="002_rasterplot_msnd2_core_"+filename+".dat"
with open(filename_rasterplot_msnd2c,"w") as results_rasterplot__msnd2c:
    for i_c in range(len(spikes_msnd2_core.t/ms)):
        results_rasterplot__msnd2c.write(str(spikes_msnd2_core.t[i_c]/ms)+","+str(spikes_msnd2_core.i[i_c])+"\n")

        
filename_rasterplot_msnd1s="003_rasterplot_msnd1_shell_"+filename+".dat"
with open(filename_rasterplot_msnd1s,"w") as results_rasterplot__msnd1s:
    for i_c in range(len(spikes_msnd1_shell.t/ms)):
        results_rasterplot__msnd1s.write(str(spikes_msnd1_shell.t[i_c]/ms)+","+str(spikes_msnd1_shell.i[i_c])+"\n")

        
filename_rasterplot_msnd2s="004_rasterplot_msnd2_shell_"+filename+".dat"
with open(filename_rasterplot_msnd2s,"w") as results_rasterplot__msnd2s:
    for i_c in range(len(spikes_msnd2_shell.t/ms)):
        results_rasterplot__msnd2s.write(str(spikes_msnd2_shell.t[i_c]/ms)+","+str(spikes_msnd2_shell.i[i_c])+"\n")        


##############################################################################
##############################################################################
##############################################################################
##############################################################################


from scipy.signal import square, sawtooth, correlate
from numpy import pi, random
import numpy
import matplotlib.pyplot as pilot


C=fr_cortex/fr_cortex[fr_cortex.argmax()]
                  
                  
V=fr_vta_dopamine/fr_vta_dopamine[fr_vta_dopamine.argmax()]
                  
T=fr_thl/fr_thl[fr_thl.argmax()]
                                  
                                  
N=fr_nacc/fr_nacc[fr_nacc.argmax()]

print("arg max kullanimi")
print(fr_cortex[fr_cortex.argmax()])







# calculate cross correlation of the two signals
xcorrCN = correlate(C, N)
xcorrVN = correlate(V, N)
xcorrTN = correlate(T, N)


recovered_phase_shift = 1-(abs(xcorrCN.argmax()-len(C)))/len(C)  

print ("Recovered phase shift: %.2f " % (recovered_phase_shift))



x1 = numpy.arange(0, len(C), 1)
x2 = numpy.arange(0, len(xcorrCN), 1)

pilot.Figure()
pilot.subplot(411)
pilot.step(x1,C)
pilot.ylabel('Korteks')
pilot.subplot(412)
pilot.step(x1,V)
pilot.ylabel('VTA')
pilot.subplot(413)
pilot.step(x1,T)
pilot.ylabel('Talamus')
pilot.subplot(414)
pilot.step(x1,N)
pilot.ylabel('NAc')
pilot.show()


pilot.Figure()
pilot.subplot(311)
pilot.step(x2,xcorrCN)
pilot.ylabel('Kortex & NAc')
pilot.subplot(312)
pilot.step(x2,xcorrVN)
pilot.ylabel('VTA & NAc')
pilot.subplot(313)
pilot.step(x2,xcorrTN)
pilot.ylabel('THL & NAc')
pilot.show()

pilot.Figure()
pilot.step(x1,fr_D1)
pilot.step(x1,fr_D2)
pilot.ylabel('D1 and D2 fr')
pilot.savefig('example2.pdf', dpi=300)
pilot.show()
#==============================================================================






#########################################
#############---------LFP-------#########
#########################################
xN_msnd1_core=5*rand(number_of_neurons_in_msnd1_core)
yN_msnd1_core=10*rand(number_of_neurons_in_msnd1_core)
zN_msnd1_core=10*rand(number_of_neurons_in_msnd1_core)
xN_msnd2_core=5*rand(number_of_neurons_in_msnd2_core)
yN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core)
zN_msnd2_core=10*rand(number_of_neurons_in_msnd2_core)
xN_msnd1_shell=5+5*rand(number_of_neurons_in_msnd1_shell)
yN_msnd1_shell=10*rand(number_of_neurons_in_msnd1_shell)
zN_msnd1_shell=10*rand(number_of_neurons_in_msnd1_shell)
xN_msnd2_shell=5+5*rand(number_of_neurons_in_msnd2_shell)
yN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell)
zN_msnd2_shell=10*rand(number_of_neurons_in_msnd2_shell)
xN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
yN_nacc_in=10*rand(number_of_neurons_in_nacc_in)
zN_nacc_in=10*rand(number_of_neurons_in_nacc_in)


array_dim=.1

xE11=5*rand()
yE11=10*rand()
zE11=10*rand()

xE21=5+5*rand()
yE21=10*rand()
zE21=10*rand()






print('Elekcrode 1 (Core): ('+str(xE11)+','+str(yE11)+','+str(zE11)+')')
print('Elekcrode 2 (Shell): ('+str(xE21)+','+str(yE21)+','+str(zE21)+')')






LFP_Is_vE1=zeros(len(trace_I_s_msnd1_core.Is[0]))*amp
LFP_Is_vE2=zeros(len(trace_I_s_msnd1_core.Is[0]))*amp

for i_lfp in range(number_of_neurons_in_msnd1_core):
    d_msnd1_core_E11=((xE11-xN_msnd1_core[i_lfp])*(xE11-xN_msnd1_core[i_lfp])+(yE11-yN_msnd1_core[i_lfp])*(yE11-yN_msnd1_core[i_lfp])+(zE11-zN_msnd1_core[i_lfp])*(zE11-zN_msnd1_core[i_lfp]))
    d_msnd2_core_E11=((xE11-xN_msnd2_core[i_lfp])*(xE11-xN_msnd2_core[i_lfp])+(yE11-yN_msnd2_core[i_lfp])*(yE11-yN_msnd2_core[i_lfp])+(zE11-zN_msnd2_core[i_lfp])*(zE11-zN_msnd2_core[i_lfp]))
    d_msnd1_shell_E11=((xE11-xN_msnd1_shell[i_lfp])*(xE11-xN_msnd1_shell[i_lfp])+(yE11-yN_msnd1_shell[i_lfp])*(yE11-yN_msnd1_shell[i_lfp])+(zE11-zN_msnd1_shell[i_lfp])*(zE11-zN_msnd1_shell[i_lfp]))
    d_msnd2_shell_E11=((xE11-xN_msnd2_shell[i_lfp])*(xE11-xN_msnd2_shell[i_lfp])+(yE11-yN_msnd2_shell[i_lfp])*(yE11-yN_msnd2_shell[i_lfp])+(zE11-zN_msnd2_shell[i_lfp])*(zE11-zN_msnd2_shell[i_lfp]))
    d_msnd1_core_E21=((xE21-xN_msnd1_core[i_lfp])*(xE21-xN_msnd1_core[i_lfp])+(yE21-yN_msnd1_core[i_lfp])*(yE21-yN_msnd1_core[i_lfp])+(zE21-zN_msnd1_core[i_lfp])*(zE21-zN_msnd1_core[i_lfp]))
    d_msnd2_core_E21=((xE21-xN_msnd2_core[i_lfp])*(xE21-xN_msnd2_core[i_lfp])+(yE21-yN_msnd2_core[i_lfp])*(yE21-yN_msnd2_core[i_lfp])+(zE21-zN_msnd2_core[i_lfp])*(zE21-zN_msnd2_core[i_lfp]))
    d_msnd1_shell_E21=((xE21-xN_msnd1_shell[i_lfp])*(xE21-xN_msnd1_shell[i_lfp])+(yE21-yN_msnd1_shell[i_lfp])*(yE21-yN_msnd1_shell[i_lfp])+(zE21-zN_msnd1_shell[i_lfp])*(zE21-zN_msnd1_shell[i_lfp]))
    d_msnd2_shell_E21=((xE21-xN_msnd2_shell[i_lfp])*(xE21-xN_msnd2_shell[i_lfp])+(yE21-yN_msnd2_shell[i_lfp])*(yE21-yN_msnd2_shell[i_lfp])+(zE21-zN_msnd2_shell[i_lfp])*(zE21-zN_msnd2_shell[i_lfp]))
    if i_lfp<50:  
        d_nacc_in_E11=((xE11-xN_nacc_in[i_lfp])*(xE11-xN_nacc_in[i_lfp])+(yE11-yN_nacc_in[i_lfp])*(yE11-yN_nacc_in[i_lfp])+(zE11-zN_nacc_in[i_lfp])*(zE11-zN_nacc_in[i_lfp]))
        d_nacc_in_E21=((xE21-xN_nacc_in[i_lfp])*(xE21-xN_nacc_in[i_lfp])+(yE21-yN_nacc_in[i_lfp])*(yE21-yN_nacc_in[i_lfp])+(zE21-zN_nacc_in[i_lfp])*(zE21-zN_nacc_in[i_lfp]))
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E11+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E11+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E11+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E11+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core[i_lfp].Is/d_msnd1_core_E21+trace_I_s_msnd2_core[i_lfp].Is/d_msnd2_core_E21+trace_I_s_msnd1_shell[i_lfp].Is/d_msnd1_shell_E21+trace_I_s_msnd2_shell[i_lfp].Is/d_msnd2_shell_E21+trace_I_s_nacc_in[i_lfp].Is/d_nacc_in_E21
    else:         
        LFP_Is_vE1=LFP_Is_vE1+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E11+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E11+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E11+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E11
        LFP_Is_vE2=LFP_Is_vE2+trace_I_s_msnd1_core.Is[i_lfp]/d_msnd1_core_E21+trace_I_s_msnd2_core.Is[i_lfp]/d_msnd2_core_E21+trace_I_s_msnd1_shell.Is[i_lfp]/d_msnd1_shell_E21+trace_I_s_msnd2_shell.Is[i_lfp]/d_msnd2_shell_E21
    




#########################################
#############------LFP (end)----#########
#########################################



#########################################
############------LFP Graphs----#########
#########################################


figure()
subplot(211)
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[0].Is / amp)
plot(trace_I_s_msnd1_core.t / ms, trace_I_s_msnd1_core[1].Is / amp)
ylabel('MSND1 Core')



subplot(212)
plot(trace_I_s_msnd1_core.t / ms, LFP_Is_vE1 /volt)
plot(trace_I_s_msnd1_core.t / ms, LFP_Is_vE2 /volt)
ylabel('LFP_Is_vE1')
xlabel('time (ms)')
tight_layout()
savefig('LFP_Is_vE1.pdf', dpi=600)



filename_LFP_Is="LFP_Is_"+filename+".dat"
with open(filename_LFP_Is,"w") as results_LFP_Is:
    for i in range(len(LFP_Is_vE1)):
 		 	results_LFP_Is.write(str(LFP_Is_vE1[i]/amp)+","+str(LFP_Is_vE2[i]/amp)+"\n")

show()



#########################################
#########--- LFP Graphs-- (end) --#######
#########################################
