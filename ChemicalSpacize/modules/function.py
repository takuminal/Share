
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



        self.mw = Descriptors.MolWt(mol)
        self.logP = Crippen.MolLogP(mol)
        self.mr =Crippen.MolMR(mol)
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




class Krevelen_sp:
    def __init__(self,smiles="C"):
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from rdkit.Chem import Crippen
        from rdkit.Chem import Draw
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        
        def pattern_count(smarts):
            return len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))

        Count  = {
        "methyl" : pattern_count("[CH3;D1]"), #メチル
        "methylene" : pattern_count("[CH2;D2]"), #メチレン
        "methyne": pattern_count("[CH;D3]"), #メチン
        "quat_carbon" : pattern_count("[C;D4]"),
        #対称性ビニルの場合にSMARTSでのカウント数を補う
        "vinyl1" : pattern_count("C=[C&D1;H2]")+pattern_count("[C&D1;H2]=[C&D1;H2]"),
        "vinyl2" : pattern_count("C=[C&D2;H1]")+pattern_count("[C&D2;H1]=[C&D2;H1]"),
        "vinyl3" : pattern_count("C=[C&D3;H0]")+pattern_count("[C&D3;H0]=[C&D3;H0]"),
        
        "cyclohexyl" : pattern_count("[C1CCCCC1]"),
        "phenyl" : pattern_count("c1ccccc1"),
        "fluoro" : pattern_count("[F]"),
        "chloro" : pattern_count("[Cl]"),
        "bromo" : pattern_count("[Br]"),
        "cyano" : pattern_count("C#N"),
        "hydroxy" : pattern_count("[C;X4][OH]"),
        "ether" : pattern_count("[C;X4]O[C;X4]"),
        "aldehyde" : pattern_count("C(=O)[H]"),
        "ketone" : pattern_count("[C]C(=O)[C]"),
        "carboxy" :pattern_count("C(=O)[OH]"),
        "ester" : pattern_count("C(=O)O[C;X4]"),
        "amino1" : pattern_count("[NH2;D1]"),
        "amino2" : pattern_count("[NH1;D2]"),
        "amino3" : pattern_count("[N;D3]"),
        }

        #原子団の寄与をdict化
        Contribute  = {
        "methyl" : [420,0,0], #メチル
        "methylene" : [270,0,0], #メチレン
        "methyne": [80,0,0], #メチン
        "quat_carbon" : [-70,0,0],
        "vinyl1" : [400,0,0],
        "vinyl2" : [200,0,0],
        "vinyl3" : [70,0,0],
        "cyclohexyl" : [1620,0,0],
        "phenyl" : [1430,110,0],
        "fluoro" : [220,0,0],
        "chloro" : [450,550,400],
        "bromo" : [550,0,0],
        "cyano" : [430,1100,2500],
        "hydroxy" : [210,500,20000],
        "ether" : [100,400,3000],
        "aldehyde" : [470,800,4500],
        "ketone" : [290,770,2000],
        "carboxy" :[530,420,10000],
        "ester" : [390,490,7000],
        "amino1" : [280,0,8400],
        "amino2" : [160,210,3100],
        "amino3" : [20,800,5000]
        }

        #分子体積の算出関数
        def mol_volume(smiles):
            from rdkit import Chem
            from rdkit.Chem import AllChem
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(mol)
            mol_vol = AllChem.ComputeMolVolume(mol)
            return mol_vol
        
              
        MolVolume = mol_volume(smiles)
        #各deltaパラメータの算出
        """Van Krevelen SPの式
        delta_d = sum(Contribute[0]*Count)/MolVolume
        delta_p = sqrt(sum((Contribute[1]**2)*Count))/MolVolume
        delta_h = sqrt(sum(Contribute[2]*Count)/MolVolume)
        sp = sqrt(delta_d**2 + delta_p**2 + delta_h**2)
        """
        ar =np.array([[Contribute[elem][0]*Count[elem],
                       ((Contribute[elem][1])**2)*Count[elem],
                       Contribute[elem][2]*Count[elem]] 
                      for elem in list(Count.keys())])
        ar_sum = ar.sum(axis = 0)
        
        delta_d = ar_sum[0]/MolVolume
        delta_p = (ar_sum[1])**(0.5)/MolVolume
        delta_h = (ar_sum[2]/MolVolume)**0.5
        delta_total = (delta_d**2 + delta_p**2 + delta_h**2)**0.5


        #アトリビュートの定義
        self.mw = Descriptors.MolWt(mol)
        self.vol = MolVolume
        self.contribute = Contribute
        self.count = Count 
        self.results =  {"delta_d":delta_d, "delta_p":delta_p, "delta_h":delta_h, "delta_total":delta_total}
    
def save_img(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(smiles)
    mol_im = Draw.MolToImage(mol=mol,size=(500,150))
    mol_im.save("images/mol_structure.png")




    
    