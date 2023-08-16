
class count_functional_groups:
    def __init__(self,smiles="C"):
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from rdkit.Chem import Crippen
        from rdkit.Chem import Draw

        mol = Chem.MolFromSmiles(smiles)
        
        def pattern_count(smarts):
            return len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))

        self.count  = {
        "carboxyl" :pattern_count("C(=O)[OH]"), # カルボキシル基のSMARTSパターン
        "carboxylate" :pattern_count("C(=O)[O-]"), # カルボキシル基のSMARTSパターン
        "hydroxyl" : pattern_count("[C;X4][OH]"),      # OH基のSMARTSパターン
        "ether" : pattern_count("[C;X4;!R]O[C;X4]"),      # 直鎖エーテル基のSMARTSパターン  
        "ester" : pattern_count("C(=O)O[C;X4]"),      # エーテル基のSMARTSパターン  
        "epoxy" : pattern_count("C1OC1"),      # エポキシ基  
        "amide" : pattern_count("C(=O)[N]"), #アミド基
        "cyano" : pattern_count("C#N"), #アミド基    
        "amine" : pattern_count("[N;X3,n]")-pattern_count("C(=O)[N]"),  # アミン基
        "ammonium" : pattern_count("[N+,n+]"),    
        "phenyl" : pattern_count("c1ccccc1"), #フェニル
        "methyl" : pattern_count("[CH3;D1]"), #メチル
        "methylene" : pattern_count("[CH2;D2]"), #メチレン
        "methyne": pattern_count("[CH;D3]"), #メチン
        "quat_carbon" : pattern_count("[C;D4]"), #4級炭素
        "vinyl" : pattern_count("[C]=[C]"), #4二重結合
        "asetylene" : pattern_count("[C]#[C]"), #三重結合
        "siloxane" :  pattern_count("O[Si](C)(C)O"),
        }

        def wiener_index(m):
            res = 0
            amat = Chem.GetDistanceMatrix(m)
            num_atoms = m.GetNumAtoms()
            for i in range(num_atoms):
                for j in range(i+1,num_atoms):
                    res += amat[i][j]
            return res

        self.mw = Descriptors.MolWt(mol)
        self.logP = Crippen.MolLogP(mol)
        self.wiener = wiener_index(mol)

        mol_im = Draw.MolToImage(mol=mol,size=(500,150))
        mol_im.save("images/mol_structure.png")


def polym(monomer):
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    from rdkit.Chem import AllChem
    
    # 反応SMARTSと反応物の定義
    rxn_smarts = "[C:1]=[C:2]>>*-[C:1][C:2]-*"
    # 反応と反応物のMoleculeオブジェクトを生成
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant = Chem.MolFromSmiles(monomer)
        products = rxn.RunReactants((reactant,))
        product_smiles = Chem.MolToSmiles(products[0][0])
    except:
        product_smiles =monomer
        print(f"reaction failed polym on {monomer}")
        pass
        
    return product_smiles


def amine_addHs(amine):
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    from rdkit.Chem import AllChem
    # 反応SMARTSと反応物の定義
    
    rxn_smarts = "[N,n:1]>>[N+,n+:1]"
    # 反応と反応物のMoleculeオブジェクトを生成
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant = Chem.MolFromSmiles(amine)
        products = rxn.RunReactants((reactant,))
        product_smiles = Chem.MolToSmiles(products[0][0])
    except:
        product_smiles = amine
        print(f"reaction failed amine_addHs on {amine}")
        pass
    return product_smiles


def acid_takeHs(acid):
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    from rdkit.Chem import AllChem

    # 反応SMARTSと反応物の定義
    rxn_smarts = "[C:1](=[O:2])[O:3]>>[C:1](=[O:2])[O-:3]"
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant = Chem.MolFromSmiles(acid)
        products = rxn.RunReactants((reactant,))
        product_smiles = Chem.MolToSmiles(products[0][0])
    except:
        product_smiles = acid
        print(f"reaction failed acid_takeHs on {acid}")
        pass
    return product_smiles


def mol_3d_viewer(smiles):
    import rdkit
    from rdkit.Chem.Draw import IPythonConsole
    from IPython.display import display, HTML
    from rdkit.Chem import Draw
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    # 3D座標を生成
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # MolオブジェクトをXYZ形式の文字列に変換する関数
    def mol_to_xyz(mol):
        num_atoms = mol.GetNumAtoms()
        xyz_string = f"{num_atoms}\nRDKit Mol to XYZ\n"
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            xyz_string += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
        return xyz_string

    # MolオブジェクトをXYZ形式に変換
    xyz_string = mol_to_xyz(mol)

    # XYZ形式の文字列をファイルに出力
    with open('output.xyz', 'w') as file:
        file.write(xyz_string)

    import py3Dmol

    v = py3Dmol.view(width=700, height=200)
    mb = Chem.MolToMolBlock(mol)
    v.addModel(mb, 'sdf')
    v.setBackgroundColor('0xeeeeee')
    v.setStyle({'stick': {}})
    v.zoomTo()
    v.show()




