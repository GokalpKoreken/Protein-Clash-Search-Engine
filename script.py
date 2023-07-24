from Bio.PDB import PDBParser
import py3Dmol 


with open("7t69.pdb", "r") as file:
    lines = file.readlines()

with open("7t69.pdb", "w") as file:
    for line in lines:
        if "ANISOU " not in line:
            file.write(line)

# Create a PDB parser
parser = PDBParser()


class Data:
    def __init__(self,atom1id,atom2id, atom1, atom2, distance, coordinates_atom1,coordinates_atom2):
        self.atom1id = atom1id
        self.atom2id = atom2id
        self.atom1 = atom1
        self.atom2 = atom2
        self.distance = distance
        self.coordinates_atom1 = coordinates_atom1
        self.coordinates_atom2 = coordinates_atom2
        


Radii = {
    "C":0.72,
    "N":0.7,
    "O":0.74,
    "H":0.34,
    "S":1
}
vdW = {
    "C":1.7,
    "N":1.55,
    "O":1.52,
    "H":1.09,
    "S":1.8
    
}




# Open the PDB file
structure = parser.get_structure("7t69", "C:/Users/gokal/OneDrive/Desktop/Protein Clash Search Engine/7t69.pdb")
atoms = structure.get_atoms()
out = "out"
i = 1
finaloutput = []


for atom1 in atoms:
    for atom2 in atoms:
        if(atom1!=atom2):       #Atomlar ayni ise atla
            cov = Radii[atom1.element] + Radii[atom2.element]
            dist = (((atom1.coord[0]-atom2.coord[0]))**2 + ((atom1.coord[1]-atom2.coord[1]))**2 + ((atom1.coord[2]-atom2.coord[2]))**2)**0.5
            if(dist<=cov or atom1.full_id[3][1] == atom2.full_id[3][1]): #Eğer covalent çapindan ufak ise veya ayni residuedalarsa atla 
                pass
            elif(atom1.full_id[3][1] == atom2.full_id[3][1] and atom1.fullid[2]!= atom2.fullid[2]): #Eğer ayni residueda ise ama zincir farklı ise hesapla
                vdw = vdW[atom1.element] + vdW[atom2.element]
                if(dist<vdw):     #Distance vdw yarıçapından den yakın ise yukarıdaki koşulları sağlamayan atom için clash var 
           
                    output = out + str(i)
                    output = Data(atom1.serial_number,atom2.serial_number,atom1.element,atom2.element, dist,atom1.coord,atom2.coord)
                    finaloutput.append(output)

            else:
                vdw = vdW[atom1.element] + vdW[atom2.element]
                if(dist<vdw):               #Distance vdw yarıçapından den yakın ise yukarıdaki koşulları sağlamayan atom için clash var
                    print(atom1.coord)  
                    output = out + str(i)
                    output = Data(atom1.serial_number,atom2.serial_number,atom1.element,atom2.element, dist,atom1.coord,atom2.coord)
                    finaloutput.append(output)


with open("7t69.pdb") as ifile:
    system = "".join([x for x in ifile])

#DİSTANCE EKLEENCEK
view = py3Dmol.view(width=700, height=700)
view.addModelsAsFrames(system)
c=1
for line in system.split("\n"):
    k=1
    split = line.split()
    if len(split) == 0 or split[0] != "ATOM":
        continue
    for data in finaloutput:
        if(data.atom1id == c or data.atom2id == c ):
            view.setStyle({'model': -1, 'serial': c+1}, {"sphere": {'color':'red','radius':0.8}})
            k=0
            break
    if(k):
        view.setStyle({'model': -1, 'serial': c+1}, {"stick": {'color':'cyan'}})
    idx = int(split[1])
    c=c+1


view.zoomTo()
view.show()