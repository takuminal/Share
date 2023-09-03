from modules import function as fc
import PySimpleGUI as sg
from rdkit import Chem
import pubchempy as pcp
import csv
import pandas as pd

sg.theme("DarkBlack1")


#官能基カウントのモジュールを空で実行し、官能基リストを作成
functional_names = fc.count_functional_groups().count.keys()

#記述子名のリスト
properties = ["Tag","ChemicalName","SMILES","MW (g/mol)","logP"]
#官能基名のリスト
functional_g = [f"{name} (mmol/g)" for name in functional_names]
#CSV出力用の２Dリスト
output_data = properties+functional_g


#過去データの取得、またはデータテーブルの初期化
try :
    #過去データの取得,DF化
    material_df = pd.read_csv("CSV/monomers.csv",encoding = "shift-jis")
    #Tagのみ
    past_tag = material_df["Tag"].tolist()
except:
    with open("CSV/monomers.csv", 'w',newline="") as file:
        writer = csv.writer(file)
        writer.writerow([output_data])
    past_tag = pd.read_csv("CSV/monomers.csv")["Tag"].tolist()


#全体のレイアウト
layout = [[sg.Text('Monomers', font=('Constantia',20,"bold"))],
          [sg.Frame("Main",[[sg.Text("化合物名 / SMILES"),
                             sg.Input(size=(40, 1), 
                                      key="-INPUT-",
                                      text_color='black',
                                      background_color='honeydew'
                                     ),
                             sg.Text("タグ"),
                             sg.Input(size=(10, 1),key = "Tag" ,
                                      text_color='black',
                                      background_color='honeydew')],
                            [sg.Radio('CAS or NAME',default=True, key="-1-", group_id='0'),
                             sg.Radio('SMILES',default=False, key="-2-", group_id='0')],
                            [sg.Text("化学反応")],
                            [sg.Checkbox("ビニル重合",key="polym", default=False)],
                            [sg.Checkbox("カルボン酸の脱プロトン化",key="Acid_DH", default=False)],
                            [sg.Checkbox("アミンのプロトン化",key="Amine_H", default=False)], 
                            [sg.Button("実行",size=(10,1)),
                             sg.Button("記録",size=(10,1)),
                             sg.Button("読込み",size=(10,1)),
                             sg.Button("削除",size=(10,1))],
                            [sg.Image(key="-IMAGE-")],

                            [sg.Text("",size=(60, 1),key = "NAME")],
                            [sg.Text("",size=(60, 1),key = "SMILES")],
                            [sg.Text("",size=(60, 1),key = "MW")],
                            [sg.Text("",size=(60, 1),key = "LogP")],
                           ]),
           
           #MainFrameの終了 / 出力テーブルの開始
           sg.Column(layout =[[sg.Text("",size=(20,1),key=name),
                              sg.Text("", key=f"{name}_")] for name in functional_names]),
           sg.Table(headings =["　Tag　"],values = [[elem] for elem in past_tag],
                    key ="output_table",
                    size=(15,30))
           ]]


# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの全画面に設定
window_size = (int(screen_width*4/5), screen_height)

# ウインドウの出現位置を指定
win_location = (0, 0)
window = sg.Window("Monomers", 
                   layout, 
                   size=window_size,
                   resizable=True,
                   location=win_location)


while True:
    event, values = window.read()
    if event == sg.WINDOW_CLOSED:
        break
    if event == "実行":
        try:
                
            if values['-1-'] == True:
    
                input_val = values["-INPUT-"]
                results = pcp.get_compounds(input_val, 'name')
                compound = results[0]
                smiles = compound.canonical_smiles
                          
            elif values['-2-'] == True:
                
                input_val = values["-INPUT-"]
                smiles = input_val
                results = pcp.get_compounds(smiles, 'smiles')
                compound = results[0]
    
            #化学反応
            
            if values["polym"] == True:
                smiles = fc.polym(smiles)
    
            if values["Acid_DH"] == True:
                smiles = fc.acid_takeHs(smiles)
    
            if values["Amine_H"] == True:
                smiles = fc.amine_addHs(smiles)                
                      
            #smilesから官能基個数のカウントモジュール実行
            fg = fc.count_functional_groups(smiles)
    
            # 画像をウィンドウ内のsg.Imageコンポーネントに更新
            window["-IMAGE-"].update("images/mol_structure.png")
            window["MW"].update(f"MW : {fg.mw} g/mol")
            window["LogP"].update(f"LogP : {fg.logP}")
            window["SMILES"].update(f"SMILES : {smiles}")
    
            for name in functional_names:
                window[name].update(f"{name} : {fg.count[name]}count")
                window[f"{name}_"].update(f"{round(1000*fg.count[name]/fg.mw,2)}mmol/g")

        except:
            pass

    if event == "記録":
        try:
                
            #追加データ行の作製
            tag = values["Tag"]
            chemical_name = values["-INPUT-"]
            new_row = [tag,chemical_name,smiles,fg.mw,fg.logP]+[1000*fg.count[name]/fg.mw for name in functional_names]
            #過去のTag名重複チェッカ－
            past_tag2 = material_df["Tag"].to_numpy()
            
            if tag in past_tag2:
                sg.popup("同じ名前が存在します",auto_close=True,auto_close_duration=1)
            else:
                new_row_df = pd.DataFrame([new_row],columns=output_data)
                material_df = pd.concat([material_df,new_row_df])
                #過去Tagデータの更新
                past_tag2 = material_df["Tag"].to_numpy()
                #Windowテーブルの更新
                window["output_table"].update([[elem] for elem in past_tag2])
                material_df.to_csv("CSV/monomers.csv",index=False,encoding ="shift-jis")
        
        except:
            pass
            
    
    if event == "削除":
        delete_Tag = values["Tag"]
        if delete_Tag in material_df["Tag"].to_list(): 
            material_df = material_df[material_df["Tag"] != delete_Tag]
            material_df = material_df[material_df["Tag"] != False ]
            past_tag2 = material_df["Tag"].to_numpy()
            window["output_table"].update([[elem] for elem in past_tag2])
            material_df.to_csv("CSV/monomers.csv",index=False,encoding ="shift-jis")

    if event == "読込み":
        Read_Name = values["Tag"]
        if Read_Name in material_df["Tag"].to_list():        
            target_row =material_df[material_df["Tag"] == Read_Name]
            smiles = target_row["SMILES"].to_numpy()[0]
            window["-INPUT-"].update(smiles)
            window["-2-"].update(True)
            
            #smilesから官能基個数のカウントモジュール実行
            fg = fc.count_functional_groups(smiles)
    
            # 画像をウィンドウ内のsg.Imageコンポーネントに更新
            window["-IMAGE-"].update("images/mol_structure.png")
            window["MW"].update(f"MW : {fg.mw} g/mol")
            window["LogP"].update(f"LogP : {fg.logP}")
            window["SMILES"].update(f"SMILES : {smiles}")
    
            for name in functional_names:
                window[name].update(f"{name} : {fg.count[name]}count")
                window[f"{name}_"].update(f"{round(1000*fg.count[name]/fg.mw,2)}mmol/g")



window.close()