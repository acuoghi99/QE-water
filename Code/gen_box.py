import numpy as np
from os import listdir, system, mkdir
import py3Dmol

m_u = 1.66e-24 # unit mass (g)
a0 = 0.53 # Bohr radius (A)

m_H = 1.008 # u
m_O = 16 # u

work_dir = './'
pack_dir = work_dir+'packmol/'

try:
    mkdir(pack_dir)
except:
    _ = None

N_H2O = 32
rho = 0.958 # g/cm^3

xyz_mol = 'water.xyz'
xyz_box = 'water_box.xyz'
pk_inp = 'water_box.inp'
pk_tol = 2.5

m_H2O = (2*m_H+m_O)*m_u
L = (N_H2O*m_H2O/rho)**(1/3)/1e-8 # Å

def gen_H2O(xyz_mol):
    r0 = [1., 1., 1.]
    theta = 104.5 #°
    d = 1 # A

    r_O = np.array(r0)
    r_H1 = r0 + np.array([d*np.cos(theta/2*np.pi/180), d*np.sin(theta/2*np.pi/180), 0])
    r_H2 = r0 + np.array([d*np.cos(theta/2*np.pi/180), -d*np.sin(theta/2*np.pi/180), 0])
    
    with open(xyz_mol, 'w') as f:
        print('3\n', file=f)
        print(f'O    {r_O[0]:.4f}    {r_O[1]:.4f}    {r_O[2]:.4f}', file=f)
        print(f'H    {r_H1[0]:.4f}    {r_H1[1]:.4f}    {r_H1[2]:.4f}', file=f)
        print(f'H    {r_H2[0]:.4f}    {r_H2[1]:.4f}    {r_H2[2]:.4f}', file=f)

if xyz_mol not in listdir(pack_dir): gen_H2O(pack_dir+xyz_mol)
        
with open(pack_dir+pk_inp, 'w') as f:
    print(f'tolerance {pk_tol:.2f}', file=f)
    print('filetype xyz', file=f)
    print(f'output {xyz_box}', file=f)
    print('', file=f)
    print(f'structure {xyz_mol}', file=f)
    print(f'  number {N_H2O}', file=f)
    print(f'inside cube 0. 0. 0. {L:.2f}', file=f)
    print('end structure', file=f)

system(f'cd {pack_dir} && packmol < {pk_inp} > packmol.out')
    
def show_system(xyz_file, L):
    p = py3Dmol.view(width=500,height=300)
    with open(xyz_file, 'r') as f:
        p.addModel(f.read(), 'xyz')
    p.setStyle({'sphere': {'radius':L/50}, 'stick':{'radius':L/100}})
    p.addLine({'start':{'x':0,'y':0,'z':0}, 'end':{'x':L,'y':0,'z':0}})
    p.addLine({'start':{'x':L,'y':0,'z':0}, 'end':{'x':L,'y':L,'z':0}})
    p.addLine({'start':{'x':L,'y':L,'z':0}, 'end':{'x':0,'y':L,'z':0}})
    p.addLine({'start':{'x':0,'y':L,'z':0}, 'end':{'x':0,'y':0,'z':0}})
    p.addLine({'start':{'x':0,'y':0,'z':0}, 'end':{'x':0,'y':0,'z':L}})
    p.addLine({'start':{'x':0,'y':0,'z':L}, 'end':{'x':L,'y':0,'z':L}})
    p.addLine({'start':{'x':L,'y':0,'z':L}, 'end':{'x':L,'y':0,'z':0}})
    p.addLine({'start':{'x':L,'y':0,'z':L}, 'end':{'x':L,'y':L,'z':L}})
    p.addLine({'start':{'x':L,'y':L,'z':L}, 'end':{'x':L,'y':L,'z':0}})
    p.addLine({'start':{'x':0,'y':0,'z':L}, 'end':{'x':0,'y':L,'z':L}})
    p.addLine({'start':{'x':0,'y':L,'z':L}, 'end':{'x':0,'y':L,'z':0}})
    p.addLine({'start':{'x':0,'y':L,'z':L}, 'end':{'x':L,'y':L,'z':L}})
    p.zoomTo()
    p.show()