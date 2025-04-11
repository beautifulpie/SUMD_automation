from Bio.PDB import PDBParser, PDBIO
import numpy as np
from Bio.PDB.Atom import Atom

input="/root/genion_test/20250318/test0317.pdb"
output="/root/genion_test/20250318/test0318_cys.pdb"

# 1. PDB 파일 파싱 (파일명은 "input.pdb"로 가정)
parser = PDBParser(QUIET=True)
structure = parser.get_structure("my_structure", input)

# 2. 모델, 체인, 잔기를 순회하면서 CYS 잔기 찾기
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() == "CYS":
                # 이미 HG가 존재하면 건너뜁니다.
                if "HG" not in residue:
                    # SG와 CA 원자가 있는지 확인
                    if "SG" in residue and "CA" in residue:
                        sg_atom = residue["SG"]
                        ca_atom = residue["CA"]
                        
                        # 3. SG와 CA 사이의 벡터 계산 및 단위 벡터 산출
                        vector = sg_atom.coord - ca_atom.coord
                        norm = np.linalg.norm(vector)
                        if norm == 0:
                            continue  # 0 나누기 방지
                        unit_vector = vector / norm
                        
                        # S-H 결합 길이 (Å): 일반적으로 약 1.33 Å
                        bond_length = 1.33
                        hg_coord = sg_atom.coord + unit_vector * bond_length
                        
                        # 4. 새로운 HG 원자 생성 (serial_number는 추후 조정 가능)
                        hg_atom = Atom(name="HG",
                                       coord=hg_coord,
                                       bfactor=sg_atom.get_bfactor(),
                                       occupancy=sg_atom.get_occupancy(),
                                       altloc=sg_atom.get_altloc(),
                                       fullname=" HG ",
                                       serial_number=0,
                                       element="H")
                        
                        # HG 원자를 잔기에 추가합니다.
                        residue.add(hg_atom)
                    else:
                        print(f"Residue {residue.get_id()}는 SG 또는 CA 원자가 없어 HG를 추가하지 않습니다.")
                        
# 5. 변경된 구조를 새로운 PDB 파일("output_with_HG.pdb")로 저장
io = PDBIO()
io.set_structure(structure)
io.save(output)
