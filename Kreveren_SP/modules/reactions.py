

def Polymerization(monomers_smiles_list):

    
    def Initiation(monomer):
        from rdkit import Chem
        from rdkit.Chem import rdChemReactions
        from rdkit.Chem import AllChem
        
        # 反応SMARTSと反応物の定義
        rxn_smarts = "[C:1]=[C:2]>>[C:1]-[C:2]-[Y]"
        # 反応と反応物のMoleculeオブジェクトを生成
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant = Chem.MolFromSmiles(monomer)
        products = rxn.RunReactants((reactant,))
        product_smiles = Chem.MolToSmiles(products[0][0])        
        return product_smiles
    
    def Propagation(active,monomer):
        from rdkit import Chem
        from rdkit.Chem import rdChemReactions
        from rdkit.Chem import AllChem
        
        # 反応SMARTSと反応物の定義
        rxn_smarts = "[C:1]-[Y].[C:2]=[C:3]>>[C:1]-[C:2]-[C:3]-[Y]"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant1 = Chem.MolFromSmiles(active)
        reactant2 = Chem.MolFromSmiles(monomer)    
        products = rxn.RunReactants((reactant1,reactant2))
        product_smiles = Chem.MolToSmiles(products[0][0])
    
        return product_smiles
    
    
    def Termination(smiles):
        from rdkit import Chem
        from rdkit.Chem import rdChemReactions
        from rdkit.Chem import AllChem
        
        # 反応SMARTSと反応物の定義
        rxn_smarts = "[C:1]-[Y]>>[C:1]"
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        reactant = Chem.MolFromSmiles(smiles)
        products = rxn.RunReactants((reactant,))
        product_smiles = Chem.MolToSmiles(products[0][0])
    
    
        return product_smiles

    import random
    feed_monomers = random.sample(monomers_smiles_list,len(monomers_smiles_list))
    initial_monomer = feed_monomers[0]
    other_monomers = feed_monomers[1:]

    radical = Initiation(initial_monomer)
    
    for monomer in other_monomers:
        radical = Propagation(radical,monomer)
    polymers = Termination(radical)

    return polymers

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


