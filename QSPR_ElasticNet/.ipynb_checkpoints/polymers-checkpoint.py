import pandas as pd 
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
import os
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib


def composition_to_function(monomer_df, composition_df):
        
    column_list = composition_df.columns.tolist()
    functional_group_df = pd.DataFrame(index = composition_df.index, 
                                       columns=list(composition_df.columns)[:2] + list(monomer_df.columns)[3:])
    
    #組成DFを1行ずつ処理
    for i in list(composition_df.index):
        #対象の1行を辞書型に変換
        composition_dict = composition_df.loc[i].to_dict()
        composition_sum = sum(list(composition_dict.values())[2:])
    
        temp=pd.DataFrame(index = [i],columns=monomer_df.columns).loc[i][3:]
        temp.fillna(0, inplace=True)
        
        #モノマー組成がある3列目以降のkeyをリスト化してループ
        for key in list(composition_dict.keys())[2:]:
            #モノマーDFのindexにkeyを照合して合致したらvalue/100を乗じてdict
            temp = monomer_df.loc[key][3:] * (composition_dict[key]/composition_sum) + temp
    
        #Seriesの情報をfunctional_group_dfに更新
        functional_group_df.loc[i]=temp
        functional_group_df.loc[i][:2]=composition_df.loc[i][:2]

    return functional_group_df





matplotlib.use('TkAgg')
monomer_df = pd.read_csv("CSV/monomers.csv",index_col = "Tag",encoding ="shift-jis")
materials = monomer_df.index.tolist()
properties = monomer_df.columns.tolist()[3:] #官能基の情報のみ


#データフレーム 
try:
    composition_df = pd.read_csv("CSV/polymer_composition.csv",encoding ="shift-jis")
    composition_df.fillna(0, inplace=True) #NaNの0フィル
except:
    composition_df = pd.DataFrame(columns=["Name","Y value"])
    composition_df.fillna(0, inplace=True) #NaNの0フィル

functional_df =  composition_to_function(monomer_df,composition_df)

#df情報をCSVに保存
functional_df.to_csv("CSV/polymer_functional_groups.csv",index=False,encoding ="shift-jis")
composition_df.to_csv("CSV/polymer_composition.csv",index=False,encoding ="shift-jis")


# GUIのレイアウト
sg.theme("DarkBlack1")

layout = [[sg.Text('Polymers', font=('Constantia',20,"bold"))],
          [sg.Text("")],
          [sg.Button("実行",size = (11,1)),
           sg.Button("データ読込",size = (11,1)),
           sg.Button("削除",size = (11,1)),
           sg.Button("グラフ描画",size = (11,1))],
          [sg.Text(' '  * 70)], #横線区切り
          [sg.Text("サンプル名",size=(15, 1)), 
           sg.Input(size=(50, 1),
                    key = "sample_name",
                    text_color='black',
                    background_color='honeydew'
                   )],
          [sg.Text("特性値",size=(15, 1)), 
           sg.Input(size=(50, 1),
                    key = "val",
                    text_color='black',
                    background_color='honeydew'
                   )],
          [sg.Text('_'  * 150)], #横線区切り
          [sg.Text('原料品名',size=(20,1)),sg.Text('仕込み量',size=(20,1))],
          [sg.Column(layout=[
              
                  [
                      sg.Combo(materials, 
                            size=(20, 1),
                            key=f'material_{n}',
                            text_color='black',
                            background_color='honeydew'
                              ),
                      sg.Input("0", 
                            size=(20, 1), 
                            key=f"feed_{n}",
                            text_color='black',
                            background_color='honeydew'
                              )
                  ] for n in range(13)],
                     
              size=(400, 400)  # 列全体のサイズ
              ),
           sg.Frame("記録データ表",
                    [[sg.Table(headings=composition_df.columns.to_list(),
                               values=composition_df.to_numpy().tolist(),
                               auto_size_columns=False,
                               col_widths=[10] * len(composition_df.columns),
                               key="output_table",
                               display_row_numbers=True,
                               vertical_scroll_only = False,
                               num_rows=min(20, 25))]],
                    vertical_alignment="top")
          ]]
# ウインドウの出現位置を指定
win_location = (0, 0)

# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの半分に設定
window_size = (screen_width,int(screen_height*0.9))

window = sg.Window("Polymers", 
                   layout, 
                   size=window_size,
                   resizable=True,
                   location=win_location)


while True:
    
    event, values = window.read()
    
    if event == sg.WINDOW_CLOSED:
        break
        
    if event == "実行":

        #出力用dict
        composition ={}

        #組成記録用dictにSample名を記録
        composition["Name"]= values["sample_name"]
        composition["Y value"]= values["val"]
        

        for i in range(13):
            material_name = values[f"material_{i}"]
            
            #materialnameに入力されている場合のみデータ参照可能
            if material_name != "":  
                #compositionに組成情報を記録
                composition[material_name]=float(values[f"feed_{i}"])
            
        #目的変数値の記録
        composition["Name"] = str(values["sample_name"])
        #目的変数値の記録
        composition["Y value"] = values["val"]

        new_name = values["sample_name"]
        #既存データの場合
        if new_name in composition_df["Name"].tolist():
            #読み取ったデータを該当の行に更新
            #行の初期化
            for elem in list(composition_df.columns)[2:]:
                composition_df.loc[composition_df["Name"] ==new_name,elem] =0
                #組成
            for key in composition.keys():
                composition_df.loc[composition_df["Name"] ==new_name,key] = composition[key]


        #新規データの場合
        else: 
            #データフレームに新しいデータ行を追加
            new_compositon_row = pd.DataFrame([composition])
            
            #Newrowをfunctional_dfに結合、インデックスは振りなおす         
            composition_df=pd.concat([composition_df,new_compositon_row],ignore_index=True)
        

        #各dfの欠損値の0フィリング
        composition_df = composition_df.fillna(0)
        functional_df = composition_to_function(monomer_df,composition_df)
        #記録リストに記録
        window["output_table"].update(composition_df.to_numpy().tolist())

        #CSVに保存
        functional_df.to_csv("CSV/polymer_functional_groups.csv",index=False,encoding ="shift-jis")
        composition_df.to_csv("CSV/polymer_composition.csv",index=False,encoding ="shift-jis")


    if event == "データ読込":

        edit_sample = values["sample_name"] #検索値をedit_sampleとして格納
        target_row = composition_df.loc[composition_df["Name"] == edit_sample] #検索値と一致する行を取得
        target_row_tr = target_row.transpose() #target_rowのdfを転置
        target_row_dict = target_row_tr.to_dict() #dictのdictとして得られる。keyはindex
        index_number = list(target_row_dict.keys())[0] #indexを数値として取得
        target_row_dict = target_row_dict[index_number] #dictの必要な部分のみ取り出して再代入
        
        window["val"].update(target_row_dict["Y value"])

        #GUIのインプットテーブルを初期化
        for i in range(13):
            window[f"material_{i}"].update("")
            window[f"feed_{i}"].update(0)
        
        #GUIのインプットテーブルの値の更新
        i = 0
        for key in list(target_row_dict.keys())[2:]:
            #値が0ではない時のみ表示
            if target_row_dict[key] !=0:
                window[f"material_{i}"].update(key)
                window[f"feed_{i}"].update(target_row_dict[key])
                i = i + 1 #次の行

    
    if event == "削除":
        try:
            delete_sample = values["sample_name"] #検索値をdelete_sampleとして格納
            composition_df = composition_df[composition_df["Name"] != delete_sample] #検索値と一致する行を削除
            functional_df = functional_df[functional_df["Name"] != delete_sample] 

            
            #削除によって空になったカラムは削除
            composition_df = composition_df.loc[:, (composition_df != 0).any(axis=0)]
            #記録リストに記録
            window["output_table"].update(composition_df.to_numpy().tolist())

            functional_df.to_csv("CSV/polymer_functional_groups.csv",index=False,encoding ="shift-jis")
            composition_df.to_csv("CSV/polymer_composition.csv",index=False,encoding ="shift-jis")
        
        except:
            pass

    if event == "グラフ描画":
            
        num_rows = 4
        num_cols = 5
        
        # ループで10個のサブプロットを生成して表示
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))
        
        for ax, elem in zip(axes.flat,functional_df.columns[2:]):
            x = functional_df[elem].tolist()
            y = functional_df["Y value"].tolist()
            ax.scatter(x, y, s=3, marker='o', color='blue')  # 仮のデータをプロット
            ax.set_title(elem,fontsize=10)
            ax.set_ylabel('Y value',fontsize=10)
        # サブプロット間のスペースを調整
        plt.tight_layout()
        plt.get_current_fig_manager().window.wm_geometry("+0+0")
        plt.show(block=False)


window.close()