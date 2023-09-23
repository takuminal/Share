
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
        "carboxy" :pattern_count("C(=O)[OH]"), # カルボキシル基のSMARTSパターン
        "carboxylate" :pattern_count("C(=O)[O-]"), # カルボキシル基のSMARTSパターン
        "hydroxy" : pattern_count("[C;X4][OH]"),      # OH基のSMARTSパターン
        "ether" : pattern_count("[C;X4;!R]O[C;X4]"),      # 直鎖エーテル基のSMARTSパターン  
        "ester" : pattern_count("C(=O)O[C;X4]"),      # エーテル基のSMARTSパターン  
        "ketone" : pattern_count("CC(=O)C"),      # エーテル基のSMARTSパターン  
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


    
def save_img(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(smiles)
    mol_im = Draw.MolToImage(mol=mol,size=(500,150))
    mol_im.save("images/mol_structure.png")





def scatter3d(x,y,z,param,tag=None,labels=["x","y","z"]):
    
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(projection='3d')

    label_x,label_y,label_z = labels
    
    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    ax.set_zlabel(label_z)

    volume = 50
    close = 0.1
    sc = ax.scatter(x,y,z,c=param, s=volume, alpha=1,cmap="coolwarm")

    ax.grid(True)
    # Add a colorbar
    colorbar = plt.colorbar(sc,shrink=0.5)
    colorbar.set_label("parameter")
    
    #データラベルの追加
    if tag is not None:
        for i in range(len(x)):
            ax.text(x[i], y[i], z[i], tag[i], fontsize=10)
        
    plt.show()

    