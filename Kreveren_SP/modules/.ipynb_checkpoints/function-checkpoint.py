
class Krevelen_sp:
    def __init__(self,smiles="C",mol_vol=None):
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
            import numpy as np 
            
            mol_vol = []
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))

            #最大5つのコンホメーションを生成し、モル体積を計算後、平均をとる。
            for i in range(3):   
                conformer_id =AllChem.EmbedMolecule(mol)
                #conformeridが0ならば、リストに追加
                if conformer_id == 0:
                    mol_vol.append(AllChem.ComputeMolVolume(mol))
                else:
                    pass
        
            result = np.mean(mol_vol)
            return result
            
        if mol_vol is None:
            mol_vol = mol_volume(smiles)
        
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
        
        delta_d = ar_sum[0]/mol_vol
        delta_p = (ar_sum[1])**(0.5)/mol_vol
        delta_h = (ar_sum[2]/mol_vol)**0.5
        delta_total = (delta_d**2 + delta_p**2 + delta_h**2)**0.5


        #アトリビュートの定義
        self.mw = Descriptors.MolWt(mol)
        self.vol = mol_vol
        self.contribute = Contribute
        self.count = Count 
        self.results =  {"delta_d":delta_d, "delta_p":delta_p, "delta_h":delta_h, "delta_total":delta_total}
    
def save_img(smiles):
    from rdkit import Chem
    from rdkit.Chem import Draw
    mol = Chem.MolFromSmiles(smiles)
    mol_im = Draw.MolToImage(mol=mol,size=(500,150))
    mol_im.save("images/mol_structure.png")


def scatter3d(x,y,z,param,tag=None):
    
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('delta_d')
    ax.set_ylabel('delta_p')
    ax.set_zlabel('delta_h')

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

def sp_distance_map(x,y,z,tag):
        import numpy as np
        
        x1,x2 = np.meshgrid(x,x.T)
        delta_d = (x1-x2)**2
        y1,y2 = np.meshgrid(y,y.T)
        delta_p = (y1-y2)**2
        z1,z2 = np.meshgrid(z,z.T)
        delta_h = (z1-z2)**2
        
        distance = np.add(delta_d,delta_p)
        distance = np.add(delta_h,distance)
        distance = distance**0.5
        
        import seaborn as sns
        import matplotlib.pyplot as plt# ヒートマップを作成
        sns.heatmap(distance, annot=False, cmap="coolwarm",xticklabels=tag, yticklabels=tag)
        # プロットを表示
        plt.show()

    