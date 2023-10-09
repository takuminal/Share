#Monomer

from modules import function as fc
import PySimpleGUI as sg
from rdkit import Chem
import pubchempy as pcp
import csv
import pandas as pd


#官能基カウントのモジュールを空で実行し、官能基リストを作成
functional_names = fc.count_functional_groups().count.keys()

#記述子名のリスト
properties = ["Tag","ChemicalName","SMILES","MW (g/mol)","logP","MR"]
#官能基名のリスト
functional_g = [name for name in functional_names]

columns_name = properties+functional_g

#過去データの取得、またはデータテーブルの初期化
try :
    #過去データの取得,DF化
    material_df = pd.read_csv("CSV/monomers.csv",encoding = "shift-jis")
    past_tag = material_df["Tag"].tolist()
except:
    material_df = pd.DataFrame(columns = columns_name)
    material_df.to_csv("CSV/monomers.csv",index=False,encoding ="shift-jis")
    past_tag = pd.read_csv("CSV/monomers.csv")["Tag"].tolist()


#全体のレイアウト
sg.theme("DarkGray6")

layout = [[sg.Text('処方管理システム / Monomer', font=("Arial",10,"bold"))],
          [sg.Frame("Main",[[sg.Text("化合物名 / SMILES"),
                             sg.Input(size=(40, 1),
                                      key="-INPUT-",
                                      background_color="gray60",text_color="darkblue"
                                     ),
                             sg.Text("タグ"),
                             sg.Input(size=(10, 1),key = "Tag",
                                      background_color="gray60",
                                      text_color="black"
                                     )
                            ],
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
                            [sg.Image("images/mol_structure.png",key="-IMAGE-")],

                            [sg.Text("SMILES",size=(10, 1)),
                             sg.Text("",size=(60, 1),key = "SMILES")],
                            [sg.Text("MW (g/mol)",size=(10, 1)),
                             sg.Text("",size=(60, 1),key = "MW")],
                            [sg.Text("LogP",size=(10, 1)),
                             sg.Text("",size=(60, 1),key = "LogP")],
                            [sg.Text("MR",size=(10, 1)),
                             sg.Text("",size=(60, 1),key = "MR")],
                           ]),
           sg.Frame("記録データ表",
                    [[sg.Table(headings=material_df.columns.to_list(),
                               values=material_df.to_numpy().tolist(),
                               auto_size_columns=False,
                               col_widths=[10] * len(material_df.columns),
                               key="output_table",
                               display_row_numbers=True,
                               vertical_scroll_only = False,
                               num_rows=min(25, 30))]],
                    vertical_alignment="top")]]

# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの全画面に設定
window_size = (screen_width, screen_height)

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
            window["MW"].update(f"{fg.mw}")
            window["LogP"].update(f"{fg.logP}")
            window["MR"].update(f"{fg.mr}")
            window["SMILES"].update(f"{smiles}")

        except:
            pass

    if event == "記録":

        #追加データ行の作製
        if values["Tag"] !="":
            tag = values["Tag"] #TagがあるときはTag名に
        else:
            tag = values["-INPUT-"] #入力がない時はINPUTを適用

        chemical_name = values["-INPUT-"]
        new_row = [tag,chemical_name,smiles,fg.mw,fg.logP,fg.mr] + [1000*fg.count[name]/fg.mw for name in functional_names]

        new_row_df = pd.DataFrame([new_row],columns = columns_name)


        #既存データの場合
        if tag in material_df["Tag"].tolist():
            #読み取ったデータを該当の行に更新
            #行の初期化
            for elem in list(material_df.columns)[1:]:
                material_df.loc[material_df["Tag"] ==tag, elem] =0
            #要素の更新
            for elem in list(material_df.columns)[1:]:
                material_df.loc[material_df["Tag"] ==tag, elem] = new_row_df.iloc[0][elem]

        #新規データの場合
        else:
            material_df = pd.concat([material_df,new_row_df])
            #Windowテーブルの更新
        window["output_table"].update(material_df.to_numpy().tolist())
        material_df.to_csv("CSV/monomers.csv",index=False,encoding ="shift-jis")

    if event == "削除":
        delete_Tag = values["Tag"]
        if delete_Tag in material_df["Tag"].to_list():
            material_df = material_df[material_df["Tag"] != delete_Tag]
            material_df = material_df[material_df["Tag"] != False ]
            past_tag2 = material_df["Tag"].to_numpy()
            window["output_table"].update(material_df.to_numpy().tolist())
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
            window["MR"].update(f"MR : {fg.mr}")
            window["SMILES"].update(f"SMILES : {smiles}")



window.close()