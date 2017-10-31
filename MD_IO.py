import numpy as np

"""
Read AMBER topology file 
"""

class AMBERFILES_TOP:

    def __init__(self, top=''):
        self.topology_file = top

        flag_atom_name = flag_charge = flag_mass =  flag_atom_type = flag_num_excl = flag_nonbonded_parm_index = flag_residue_label = flag_residue_pointer = flag_bond_const = \
        flag_bond_value = flag_angle_const = flag_angle_value = flag_dihedral_const = flag_dihedral_peri = \
        flag_dihedral_phase = flag_LJ_A = flag_LJ_B = flag_bond_inc_h = flag_bond_wo_h = \
        flag_angle_inc_h = flag_angle_wo_h = flag_dihedral_inc_h = flag_dihedral_wo_h = count = flag_excl_list = flag_hbond_A = flag_hbond_B = flag_atom_type_name = 0

        charge_list = mass_list = bond_const_list = bond_value_list = angle_const_list = angle_value_list = dihedral_const_list = np.array([], dtype=float) 
        dihedral_phase_list = LJ_A_list = LJ_B_list = coor_list = excl_list = hbond_A_list = hbond_B_list = np.array([], dtype=float) 

        atom_type_list = num_excl_list = residue_pointer = nonbonded_parm_index_list = dihedral_peri_list = np.array([], dtype=int) 
        bond_inc_h_list = bond_wo_h_list = angle_inc_h_list = angle_wo_h_list = dihedral_inc_h_list = dihedral_wo_h_list = np.array([], dtype=int)
    
        atom_name = residue_label = atom_type_name = np.array([], dtype = str)

        box = angle = np.zeros(3)

        for line in open(self.topology_file):
            line = line.split()
            if(len(line) >=1):
                if(line [0] == '%FLAG' and line[1] == 'ATOM_NAME'):
                    flag_atom_name = 1
                elif(line [0] == '%FLAG' and line[1] == 'CHARGE'):
                    flag_atom_name = 0
                    flag_charge = 1
                elif(line [0] == '%FLAG' and line[1] == 'ATOMIC_NUMBER'):
                    flag_charge = 0
                elif(line [0] == '%FLAG' and line[1] == 'MASS'):
                    flag_mass = 1
                elif(line [0] == '%FLAG' and line[1] == 'ATOM_TYPE_INDEX'):
                    flag_mass = 0
                    flag_atom_type = 1
                elif(line [0] == '%FLAG' and line[1] == 'NUMBER_EXCLUDED_ATOMS'):
                    flag_atom_type = 0
                    flag_num_excl = 1
                elif(line [0] == '%FLAG' and line[1] == 'NONBONDED_PARM_INDEX'):
                    flag_num_excl = 0
                    flag_nonbonded_parm_index = 1
                elif(line [0] == '%FLAG' and line[1] == 'RESIDUE_LABEL'):
                    flag_residue_label = 1
                    flag_nonbonded_parm_index = 0
		elif(line [0] == '%FLAG' and line[1] == 'RESIDUE_POINTER'):
                    flag_residue_label = 0
                    flag_residue_pointer = 1
                elif(line [0] == '%FLAG' and line[1] == 'BOND_FORCE_CONSTANT'):
                    flag_residue_pointer = 0
                    flag_bond_const = 1
                elif(line [0] == '%FLAG' and line[1] == 'BOND_EQUIL_VALUE'):
                    flag_bond_const = 0
                    flag_bond_value = 1
                elif(line [0] == '%FLAG' and line[1] == 'ANGLE_FORCE_CONSTANT'):
                    flag_angle_const = 1
                    flag_bond_value = 0
                elif(line [0] == '%FLAG' and line[1] == 'ANGLE_EQUIL_VALUE'):
                    flag_angle_const = 0
                    flag_angle_value = 1
                elif(line [0] == '%FLAG' and line[1] == 'DIHEDRAL_FORCE_CONSTANT'):
                    flag_dihedral_const = 1
                    flag_angle_value = 0
                elif(line [0] == '%FLAG' and line[1] == 'DIHEDRAL_PERIODICITY'):
                    flag_dihedral_const = 0
                    flag_dihedral_peri = 1
                elif(line [0] == '%FLAG' and line[1] == 'DIHEDRAL_PHASE'):
                    flag_dihedral_phase = 1
                    flag_dihedral_peri = 0
                elif(line [0] == '%FLAG' and line[1] == 'SCEE_SCALE_FACTOR'):
                    flag_dihedral_phase = 0
                elif(line [0] == '%FLAG' and line[1] == 'LENNARD_JONES_ACOEF'):
                    flag_LJ_A = 1
                elif(line [0] == '%FLAG' and line[1] == 'LENNARD_JONES_BCOEF'):
                    flag_LJ_B = 1
                    flag_LJ_A = 0
                elif(line [0] == '%FLAG' and line[1] == 'BONDS_INC_HYDROGEN'):
                    flag_LJ_B = 0
                    flag_bond_inc_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'BONDS_WITHOUT_HYDROGEN'):
                    flag_bond_inc_h = 0
                    flag_bond_wo_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'ANGLES_INC_HYDROGEN'):
                    flag_bond_wo_h = 0
                    flag_angle_inc_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'ANGLES_WITHOUT_HYDROGEN'):
                    flag_angle_inc_h = 0
                    flag_angle_wo_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'DIHEDRALS_INC_HYDROGEN'):
                    flag_angle_wo_h = 0
                    flag_dihedral_inc_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'DIHEDRALS_WITHOUT_HYDROGEN'):
                    flag_dihedral_inc_h = 0
                    flag_dihedral_wo_h = 1
                elif(line [0] == '%FLAG' and line[1] == 'EXCLUDED_ATOMS_LIST'):
                    flag_dihedral_wo_h = 0
                    flag_excl_list = 1
                elif(line [0] == '%FLAG' and line[1] == 'HBOND_ACOEF'):
                    flag_excl_list = 0
                    flag_hbond_A = 1
                elif(line [0] == '%FLAG' and line[1] == 'HBOND_BCOEF'):
                    flag_hbond_A = 0
                    flag_hbond_B = 1
                elif(line [0] == '%FLAG' and line[1] == 'HBCUT'):
                    flag_hbond_B = 0
                elif(line [0] == '%FLAG' and line[1] == 'AMBER_ATOM_TYPE'):
                    flag_atom_type_name = 1
                elif(line [0] == '%FLAG' and line[1] == 'TREE_CHAIN_CLASSIFICATION'):
                    flag_atom_type_name = 0
                else:
                    a=0

                if(flag_atom_name == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = str(line[i])
                    atom_name = np.append(atom_name, line)
                elif(flag_charge == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    charge_list = np.append(charge_list, line)
                elif(flag_mass == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    mass_list = np.append(mass_list, line)
                elif(flag_atom_type == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    atom_type_list = np.append(atom_type_list, line)
                elif(flag_num_excl == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    num_excl_list = np.append(num_excl_list, line)
                elif(flag_nonbonded_parm_index == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    nonbonded_parm_index_list = np.append(nonbonded_parm_index_list, line)
                elif(flag_residue_label == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = str(line[i])
                    residue_label = np.append(residue_label, line)
                elif(flag_residue_pointer == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    residue_pointer = np.append(residue_pointer, line)
                elif(flag_bond_const == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    bond_const_list = np.append(bond_const_list, line)
                elif(flag_bond_value == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    bond_value_list = np.append(bond_value_list, line)
                elif(flag_angle_const == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    angle_const_list = np.append(angle_const_list, line)
                elif(flag_angle_value == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    angle_value_list = np.append(angle_value_list, line)
                elif(flag_dihedral_const == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    dihedral_const_list = np.append(dihedral_const_list, line)
                elif(flag_dihedral_peri == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    dihedral_peri_list = np.append(dihedral_peri_list, line)
                elif(flag_dihedral_phase == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    dihedral_phase_list = np.append(dihedral_phase_list, line)
                elif(flag_LJ_A == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    LJ_A_list = np.append(LJ_A_list, line)
                elif(flag_LJ_B == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                    for i in range(len(line)):
                        line[i] = float(line[i])
                    LJ_B_list = np.append(LJ_B_list, line)
                elif(flag_bond_inc_h == 1 and len(line) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    bond_inc_h_list = np.append(bond_inc_h_list, line)
                elif(flag_bond_wo_h == 1 and len(line) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    bond_wo_h_list = np.append(bond_wo_h_list, line)
                elif(flag_angle_inc_h == 1 and len(line) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    angle_inc_h_list = np.append(angle_inc_h_list, line)
                elif(flag_angle_wo_h == 1 and len(line) >= 1  and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    angle_wo_h_list = np.append(angle_wo_h_list, line)
                elif(flag_dihedral_inc_h == 1 and len(line) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    dihedral_inc_h_list = np.append(dihedral_inc_h_list, line)
                elif(flag_dihedral_wo_h == 1 and len(line) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    dihedral_wo_h_list = np.append(dihedral_wo_h_list, line)
                elif(flag_excl_list == 1 and len(line) >=1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = int(line[i])
                    excl_list = np.append(excl_list, line)
                elif(flag_hbond_A == 1 and len(line) >=1):
                    if(flag_hbond_A == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                        for i in range(len(line)):
                            line[i] = int(line[i])
                        hbond_A_list = np.append(hbond_A_list, line)
                elif(flag_hbond_B == 1 and len(line) >=1):
                    if(flag_hbond_B == 1 and len(line[0]) >= 14 and not 'FORMAT' in line[0]):
                        for i in range(len(line)):
                            line[i] = int(line[i])
                        hbond_B_list = np.append(hbond_B_list, line)
                elif(flag_atom_type_name == 1 and len(line[0]) >= 1 and not 'FORMAT' in line[0] and not 'FLAG' in line[0]):
                    for i in range(len(line)):
                        line[i] = str(line[i])
                    atom_type_name = np.append(atom_type_name, line)

        self._atom_name = atom_name
	self._charge_list = charge_list
        self._mass_list = mass_list
        self._atom_type_list = atom_type_list
        self._num_excl_list = num_excl_list
        self._residue_label = residue_label
        self._residue_pointer = residue_pointer
        self._nonbonded_parm_index_list = nonbonded_parm_index_list
        self._bond_const_list = bond_const_list
        self._bond_value_list = bond_value_list
        self._angle_const_list = angle_const_list
        self._angle_value_list = angle_value_list
        self._dihedral_const_list = dihedral_const_list
        self._dihedral_peri_list = dihedral_peri_list
        self._dihedral_phase_list = dihedral_phase_list
        self._LJ_A_list = LJ_A_list
        self._LJ_B_list = LJ_B_list
        self._bond_inc_h_list = bond_inc_h_list
        self._bond_wo_h_list =  bond_wo_h_list
        self._angle_inc_h_list = angle_inc_h_list
        self._angle_wo_h_list = angle_wo_h_list
	self._dihedral_inc_h_list = dihedral_inc_h_list
        self._dihedral_wo_h_list = dihedral_wo_h_list
	self._excl_list = excl_list
        self._hbond_A_list = hbond_A_list
        self._hbond_B_list = hbond_B_list
        self._atom_type_name = atom_type_name

class resinfo(AMBERFILES_TOP):
	def resid(self, resid):
            self.resid = resid

            return self.resid

        def resname(self, resid):
            self.resname = self._residue_label[resid-1]

            return self.resname

	def atoms(self, resid):
            self.atoms = self._atom_name[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]

            return self.atoms

        def atoms_type(self, resid):
            self.atoms_type = self._atom_type_name[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]

            return self.atoms_type

        def charges(self, resid):
            self.charges = self._charge_list[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]/18.2223

            return self.charges

        def bondterms(self, resid):
            begin = self._residue_pointer[resid-1]
            end = self._residue_pointer[resid]

            self.bondterms = np.array([])

            for i in range(0, len(self._bond_inc_h_list), 3):
                atom_a_id = self._bond_inc_h_list[i] / 3 + 1
                atom_b_id = self._bond_inc_h_list[i+1] / 3 + 1
                bond_id = self._bond_inc_h_list[i+2]
                if atom_a_id >= begin and atom_a_id < end and atom_b_id >= begin:
                    self.bondterms = np.append(self.bondterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) + ') ' + str(self._bond_const_list[bond_id-1]) + " " + str(self._bond_value_list[bond_id-1]))

            for i in range(0, len(self._bond_wo_h_list), 3):
                atom_a_id = self._bond_wo_h_list[i] / 3 + 1
                atom_b_id = self._bond_wo_h_list[i+1] / 3 + 1
                bond_id = self._bond_wo_h_list[i+2]
                if atom_a_id >= begin and atom_a_id < end and atom_b_id >= begin:
                    self.bondterms = np.append(self.bondterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) +') ' + str(self._bond_const_list[bond_id-1]) + " " + str(self._bond_value_list[bond_id-1]))

            return self.bondterms

        def angleterms(self, resid):
            begin = self._residue_pointer[resid-1]
            end = self._residue_pointer[resid]
            self.angleterms = np.array([])

            for i in range(0, len(self._angle_inc_h_list), 4):
                atom_a_id = self._angle_inc_h_list[i] / 3 + 1
                atom_b_id = self._angle_inc_h_list[i+1] / 3 + 1
                atom_c_id = self._angle_inc_h_list[i+2] / 3 + 1
                angle_id = self._angle_inc_h_list[i+3]
                if atom_a_id >= begin and atom_a_id < end:
                    self.angleterms = np.append(self.angleterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) + ')-' + str(self._atom_name[atom_c_id-1]) + '(' + str(atom_c_id) + ') ' + str(self._angle_const_list[angle_id-1]) + " " + str(self._angle_value_list[angle_id-1]))

            for i in range(0, len(self._angle_wo_h_list), 4):
                atom_a_id = self._angle_wo_h_list[i] / 3 + 1
                atom_b_id = self._angle_wo_h_list[i+1] / 3 + 1
                atom_c_id = self._angle_wo_h_list[i+2] / 3 + 1
                angle_id = self._angle_wo_h_list[i+3]
                if atom_a_id >= begin and atom_a_id < end:
                    self.angleterms = np.append(self.angleterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) + ')-' + str(self._atom_name[atom_c_id-1]) + '(' + str(atom_c_id) + ') ' + str(self._angle_const_list[angle_id-1]) + " " + str(self._angle_value_list[angle_id-1]))

            return self.angleterms

        def dihedralterms(self, resid):
            begin = self._residue_pointer[resid-1]
            end = self._residue_pointer[resid]
            self.dihedralterms = np.array([])
            for i in range(0, len(self._dihedral_inc_h_list), 5):
                atom_a_id = self._dihedral_inc_h_list[i] / 3 + 1
                atom_b_id = self._dihedral_inc_h_list[i+1] / 3 + 1
                atom_c_id = self._dihedral_inc_h_list[i+2] / 3 + 1
                atom_d_id = self._dihedral_inc_h_list[i+3] / 3 + 1
                dihedral_id = self._dihedral_inc_h_list[i+4]
                if atom_a_id >= begin and atom_a_id < end:
                    self.dihedralterms = np.append(self.dihedralterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) + ')-' + str(self._atom_name[atom_c_id-1]) + '(' + str(atom_c_id) + ')-' + str(self._atom_name[atom_d_id-1]) + '(' + str(atom_d_id) + ') ' + str(self._dihedral_const_list[dihedral_id-1]) + " " + str(self._dihedral_phase_list[dihedral_id-1]) + " " + str(self._dihedral_peri_list[dihedral_id-1]))
   
            for i in range(0, len(self._dihedral_wo_h_list), 5):
                atom_a_id = self._dihedral_wo_h_list[i] / 3 + 1
                atom_b_id = self._dihedral_wo_h_list[i+1] / 3 + 1
                atom_c_id = self._dihedral_wo_h_list[i+2] / 3 + 1
                atom_d_id = self._dihedral_wo_h_list[i+3] / 3 + 1
                dihedral_id = self._dihedral_wo_h_list[i+4]
                if atom_a_id >= begin and atom_a_id < end:
                    self.dihedralterms = np.append(self.dihedralterms, str(self._atom_name[atom_a_id-1]) + '(' + str(atom_a_id) + ')-' + str(self._atom_name[atom_b_id-1]) + '(' + str(atom_b_id) + ')-' + str(self._atom_name[atom_c_id-1]) + '(' + str(atom_c_id) + ')-' + str(self._atom_name[atom_d_id-1]) + '(' + str(atom_d_id) + ') ' + str(self._dihedral_const_list[dihedral_id-1]) + " " + str(self._dihedral_phase_list[dihedral_id-1]) + " " + str(self._dihedral_peri_list[dihedral_id-1]))            
    
            return self.dihedralterms

        def summary(self, resid):
            print "==================================="
            print "Summary of residue " + str(resid)
	    print "-----------------------------------"
            print "Residue Name: " + str(self._residue_label[resid-1])
            print "Net Charge: " + str(sum(self._charge_list[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]/18.2223))
            print "Number of Atom(s): " + str(len(self._atom_name[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]))
            print "Atom(s): " + str(self._atom_name[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1])
            print "Atom type: " + str(self._atom_type_name[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1])
            print "Partial Charges: " + str(self._charge_list[self._residue_pointer[resid-1]-1:self._residue_pointer[resid]-1]/18.2223)


