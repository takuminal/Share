#Distribute

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
import os
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns


matplotlib.use('TkAgg')

polymer_df = pd.read_csv("CSV/polymer_functional_groups.csv",index_col="Name",encoding ="shift-jis")#ポリマーデータの読み込み
additives_df = pd.read_csv("CSV/additives.csv",index_col="Name",encoding ="shift-jis")
results_df = pd.read_csv("CSV/results.csv",index_col="Name",encoding ="shift-jis")
test_name = results_df.index.tolist()
solvent_names = ["PW","MIBK","Toluene","MEK","Butyl Acetate","iBuOH","BuOH","IPA"]



try:
    brend_df = pd.read_csv("CSV/brend.csv",encoding ="shift-jis")
    log_df = pd.read_csv("CSV/brend_log.csv",encoding ="shift-jis")
except:
    col_1 = ["Name"]
    col_2 = polymer_df.columns[3:].tolist()
    col_3 = additives_df.columns[2:].tolist()
    col_4 = results_df.index.tolist()

    brend_df = pd.DataFrame(columns=col_1+col_2+col_3+col_4)
    log_df = pd.DataFrame([])


polymer_names = polymer_df.index.tolist()
additive_names = additives_df.index.tolist()

#処方出力用のデータフレーム
shoho_head = [f"{i}" for i in range(10)]
shoho_df=pd.DataFrame([])

# GUIのレイアウト
sg.theme("DefaultNoMoreNagging")

layout = [[sg.Text("処方管理システム / 配合", font=('System',10,"bold"))],
          [sg.Frame("Main",
                    layout = [
                        [sg.Button("計算",size = (10,1)),
                         sg.Button("削除",size = (10,1)),
                         sg.Button("読み込み",size = (10,1)),
                         sg.Button("登録",size = (10,1)),

                        ],

                        [sg.Text('_'  * 70)], #横線区切り
                        [sg.Text("サンプル名",size=(15, 1)),
                         sg.Combo(brend_df["Name"].tolist(),
                                  size=(50, 1),
                                  key = "sample_name")],

                        [sg.Text("濃度",size=(15, 1)),
                         sg.Input(size=(10, 1),
                                  default_text="10",
                                  key = "conc")],
                        [sg.Text("スケール",size=(15, 1)),
                         sg.Input(size=(10, 1),
                                  default_text="10",
                                  key = "scale")],
                        [sg.Text('_'  * 70)], #横線区切り
                        [sg.Text('主剤',size=(20,1)),
                         sg.Text('固形%',size=(10,1)),
                         sg.Text('固形量',size=(10,1),),
                         sg.Text('仕込み量',size=(10,1)),

                        ],
                        [sg.Column(layout=[
                            [sg.Combo(polymer_names, size=(20, 1),key=f'main_{n}'),
                             sg.Input("0", size=(10, 1),key=f"main_solid_ratio_{n}"),
                             sg.Text("0", size=(10, 1),key=f"main_solid_{n}"),
                             sg.Input("0", size=(10, 1),key=f"main_feed_{n}"
                                      ,background_color="gray90"),
                            ] for n in range(2)],size=(400, 50)  # 列全体のサイズ
                                        )
                        ],
                        [sg.Text('添加剤',size=(20,1))],
                        [sg.Column(layout=[
                            [sg.Combo(additive_names, size=(20, 1),key=f'ad_{n}'),
                             sg.Input("0", size=(10, 1), key=f"ad_solid_ratio_{n}"),
                             sg.Text("0", size=(10, 1), key=f"ad_solid_{n}"),
                             sg.Input("0", size=(10, 1), key=f"ad_feed_{n}"
                                      ,background_color="gray90"),


                            ] for n in range(5)],
                                   size=(400,100)  # 列全体のサイズ
                                  )
                              ],
                        [sg.Text("希釈溶媒",size=(20,1))],
                        [sg.Combo(solvent_names,size=(20,1),key="solvent"),
                         sg.Text("",size=(20,1)),
                         sg.Input("",size=(10,1),key="solvent_wt",background_color="gray90")],

                        [sg.Text('試験名',size=(20,1)),sg.Text('値',size=(20,1))],
                        [sg.Column(layout=[
                            [sg.Combo(test_name,size=(20, 1),key=f'res_name_{n}'),
                             sg.Input("0", size=(10, 1), key=f"res_value_{n}")
                            ] for n in range(3)],
                                   size=(400, 100)  # 列全体のサイズ
                                  )
                        ]
                    ]),
          sg.Frame("記録データ表",
                    [[sg.Table(headings=brend_df.columns.to_list(),
                      values=brend_df.to_numpy().tolist(),
                      auto_size_columns=False,
                      col_widths=[10] * len(brend_df.columns),
                      key="output_table",
                      display_row_numbers=True,
                      vertical_scroll_only = False,
                      num_rows=min(20, 15))],

                     [sg.Text("処方名",size=(15, 1)),
                      sg.Input("",size=(30, 1),key = "shoho_name")],
                     [sg.Text("サンプル選択",size=(15, 1)),
                      sg.Combo(brend_df["Name"].tolist(),size=(30, 1),key = "shoho-ad"),
                      sg.Button("処方に追加")],

                     [sg.Table(headings=shoho_head,
                      values=shoho_df.to_numpy().tolist(),
                      auto_size_columns=False,
                      col_widths=[10] * len(shoho_df.columns),
                      key="output_table2",
                      display_row_numbers=True,
                      vertical_scroll_only = False,
                      num_rows=min(20, 15))]
                    ],
                   vertical_alignment="top")
          ]
         ]



# ウインドウの出現位置を指定
win_location = (0, 0)

# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの半分に設定
window_size = (screen_width,int(screen_height*0.9))

window = sg.Window("Brend",
                   layout,
                   size=window_size,
                   resizable=True,
                   location=win_location)
while True:
    event, values = window.read()
    if event == sg.WINDOW_CLOSED:
        break

    if event == "計算":

        conc = float(values["conc"])
        scale = float(values["scale"])
        total_weight = scale
        total_solid  = scale*conc/100
        wt_except_solvent =0

        function_1=pd.DataFrame({key:[0] for key in polymer_df.columns[3:]})
        for n in range(2):
            name = values[f'main_{n}']
            if name !="":

                solid_ratio = float(values[f"main_solid_ratio_{n}"])
                solid = (solid_ratio/100)*total_solid
                nv =polymer_df["Conc"][name]
                feed = round(solid/(nv/100),3)
                wt_except_solvent = wt_except_solvent + feed
                function_1 = polymer_df.loc[name][3:]*solid_ratio + function_1

            elif name == "":
                solid = 0
                feed = 0

            window[f"main_solid_{n}"].update(f"{solid} g")
            window[f"main_feed_{n}"].update(f"{feed} g")

        function_2=pd.DataFrame({key:[0] for key in additives_df.columns[2:]})
        for n in range(5):
            name = values[f'ad_{n}']
            if name !="":

                solid_ratio = float(values[f"ad_solid_ratio_{n}"])
                solid = (solid_ratio/100)*total_solid
                nv =additives_df["Conc"][name]
                feed = round(solid/(nv/100),3)
                wt_except_solvent = round(wt_except_solvent + feed,3)
                function_2 = additives_df.loc[name][2:]*solid_ratio + function_2

            elif name == "":
                solid = 0
                feed = 0

            window[f"ad_feed_{n}"].update(f"{feed} g")
            window[f"ad_solid_{n}"].update(f"{solid} g")

        solvent_wt = total_weight - wt_except_solvent
        window["solvent_wt"].update(f"{solvent_wt} g")

        results = {} #試験結果の格納場所
        for n in range(3):
            target_name = values[f'res_name_{n}']
            target_value= values[f"res_value_{n}"]
            if target_name!="":
                results[target_name]=[target_value]

        results_df = pd.DataFrame(results)

        #------------------------------------------------------------------
        new_row_name=pd.DataFrame({"Name":[values["sample_name"]]})
        #名前、主剤、添加剤、結果をconcatして一行にする
        new_row = pd.concat([new_row_name,function_1,function_2,results_df],axis=1)




    if event == "登録":

        #現在のUIに表示されている情報を読み込み、new_logとして新しいlogの行を作成。
        out = values
        out2={}
        for key,val in out.items():
            out2[key]=[val]
        new_log = pd.DataFrame(out2)

        #データの記録-------------------------------------------------------
        target_name = values["sample_name"]

        #新規データ
        if target_name not in brend_df["Name"].to_list():
            #brend_dfに新しい行を追加
            brend_df = pd.concat([brend_df,new_row])
            brend_df = brend_df.sort_values("Name")
            brend_df = brend_df.sort_values("Name")
            window["output_table"].update(brend_df.to_numpy().tolist())
            brend_df.to_csv("CSV/brend.csv",index=False,encoding ="shift-jis")

            log_df = pd.concat([new_log,log_df]).fillna("")
            log_df = log_df.sort_values("sample_name")
            log_df = log_df.reset_index(drop=True)
            log_df.to_csv("CSV/brend_log.csv",index=False,encoding ="shift-jis")

            #-------------------------------------------------------------------------
        #既存データ
        elif target_name in brend_df["Name"].to_list():

            target_name = values["sample_name"]
            brend_df = brend_df[brend_df["Name"]!=target_name]
            log_df = log_df[log_df["sample_name"]!=target_name]

            brend_df = pd.concat([brend_df,new_row])
            brend_df = brend_df.sort_values("Name")
            brend_df = brend_df.reset_index(drop=True)
            brend_df.to_csv("CSV/brend.csv",index=False,encoding ="shift-jis")

            log_df = pd.concat([new_log,log_df]).fillna("")
            log_df = log_df.sort_values("sample_name")
            log_df = log_df.reset_index(drop=True)
            log_df.to_csv("CSV/brend_log.csv",index=False,encoding ="shift-jis")

            window["output_table"].update(brend_df.to_numpy().tolist())

        else:
            pass

        window["shoho-ad"].update(values = log_df["sample_name"].tolist())


    if event =="読み込み":
        target_name = values["sample_name"]
        target_row = log_df[log_df["sample_name"]==target_name]
        target_row=target_row.set_index("sample_name")
        target_row.pop("output_table")
        target_row.pop("output_table2")
        target_row = target_row.fillna("")

        for col in target_row.columns:
            val = target_row[col][target_name]
            window[col].update(val)

    if event == "削除":
        del_name = values["sample_name"]

        brend_df = brend_df[brend_df["Name"]!=del_name]
        brend_df = brend_df.sort_values("Name")


        log_df = log_df[log_df["sample_name"]!=del_name]
        log_df = log_df.sort_values("sample_name")

        brend_df.to_csv("CSV/brend.csv",index=False,encoding ="shift-jis")
        log_df.to_csv("CSV/brend_log.csv",index=False,encoding ="shift-jis")

        window["output_table"].update(brend_df.to_numpy().tolist())


    if event == "処方に追加":

        #log_dfのターゲット行を辞書型に戻す
        target_name = values["shoho-ad"]
        target_row = log_df[log_df["sample_name"] == target_name].index
        target_row_index = target_row.tolist()[0]
        target_dict = log_df.loc[target_row_index].to_dict()

        target_dict_re={}
        for key,val in target_dict.items():
            target_dict_re[key]=[val]
        new_log = pd.DataFrame(target_dict_re)
        exp_sheet = pd.DataFrame([])
        exp_sheet = pd.concat([exp_sheet,new_log.T],axis = 1)

        shoho={}

        #名前情報の付加
        shoho["サンプル名"]= [exp_sheet.loc["sample_name"][0],"",""]
        shoho["濃度 (%)"]= [exp_sheet.loc["conc"][0],"",""]
        shoho["スケール (g)"]= [exp_sheet.loc["scale"][0],"",""]

        #主剤
        for i in range(2):
            col_o = f"主剤{i}"
            col_a = exp_sheet.loc[f"main_{i}"][0]
            col_b = exp_sheet.loc[f"main_feed_{i}"][0]
            shoho[col_o] = [col_a,col_b,""]

        #添加剤
        for i in range(5):

            col_o =f"添加剤{i}"
            col_a = exp_sheet.loc[f"ad_{i}"][0]
            col_b = exp_sheet.loc[f"ad_feed_{i}"][0]
            shoho[col_o] = [col_a,col_b,""]

        #希釈溶剤
        col_o = "溶剤"
        col_a = exp_sheet.loc["solvent"][0]
        col_b = exp_sheet.loc["solvent_wt"][0]
        shoho[col_o] = [col_a,col_b,""]

        shoho_name = values["shoho_name"]

        new_shoho_col = pd.DataFrame(shoho)
        shoho_df = pd.concat([shoho_df,new_shoho_col],axis = 0)
        shoho_df.to_csv("処方/"+shoho_name+".csv",encoding ="shift-jis")
        window["output_table2"].update(shoho_df.to_numpy().tolist())

window.close()
