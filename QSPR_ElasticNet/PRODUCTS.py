import pandas as pd 
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
import os
import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

df = pd.read_csv("CSV/materials.csv",index_col = "Tag")
materials = df.index.tolist()
properties = df.columns.tolist()[3:] #官能基の情報のみ


#機械学習用データフレーム 
try:
    functional_df = pd.read_csv("CSV/functional_groups.csv")
    composition_df = pd.read_csv("CSV/composition.csv")
    

except:
    funcional_df = pd.DataFrame(columns=["Name","Object Val"]+properties)
    composition_df = pd.DataFrame(columns=["Name","Object Val"])

    

# GUIのレイアウト
sg.theme("DarkBlack1")

layout = [[sg.Text('PRODUCTS', font=('Constantia',20,"bold"))],
          [sg.Text("")],
          [sg.Button("実行",size = (11,1)),
           sg.Button("データ読込",size = (11,1)),
           sg.Button("削除",size = (11,1)),
           sg.Button("グラフ描画",size = (11,1)),
           sg.Button("保存",size = (11,1))],
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
          [sg.Text('_'  * 70)], #横線区切り
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

window = sg.Window("PRODUCTS", 
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
        output={} 
        composition ={}

        #組成記録用dictにSample名を記録
        composition["Name"]= values["sample_name"]
        composition["Object Val"]= values["val"]
        
        
        for property in properties:
            calc_val = []
            feed_sum = []
            
            for i in range(13):
                material_name = values[f"material_{i}"]
                if material_name != "":  #materialnameに入力されている場合のみデータ参照可能
                    feed = float(values[f"feed_{i}"])
                    calc_val.append(df.loc[material_name,property]*feed)
                    feed_sum.append(feed)

                    #compositionに組成情報を記録
                    composition[material_name]=values[f"feed_{i}"]

                
                #目的変数値の記録
                output["Name"] = str(values["sample_name"])
                #目的変数値の記録
                output["Object Val"] = values["val"]
                #官能基情報の記録
                output[property] = sum(calc_val)/sum(feed_sum)

       
        new_name = values["sample_name"]
        #既存データの場合
        if (new_name in functional_df["Name"].tolist()) and (new_name in composition_df["Name"].tolist()):
            #読み取ったデータを該当の行に更新
                #官能基
            for key in output.keys():
                functional_df.loc[functional_df["Name"] == new_name,key] = output[key]
                
                #組成
            for key in composition.keys():
                composition_df.loc[composition_df["Name"] ==new_name,key] = composition[key]

        #新規データの場合
        else: 
            #データフレームに新しいデータ行を追加
            new_row = pd.DataFrame([output])
            new_compositon_row = pd.DataFrame([composition])
            
            #Newrowをfunctional_dfに結合、インデックスは振りなおす
            functional_df = pd.concat([functional_df,new_row],ignore_index=True)            
            composition_df=pd.concat([composition_df,new_compositon_row])

        #各dfの欠損値の0フィリング
        functional_df = functional_df.fillna(0)
        composition_df = composition_df.fillna(0)

        #記録リストに記録
        window["output_table"].update(composition_df.to_numpy().tolist())


    if event == "データ読込":

        edit_sample = values["sample_name"] #検索値をedit_sampleとして格納
        target_row = composition_df.loc[composition_df["Name"] == edit_sample] #検索値と一致する行を取得
        target_row_tr = target_row.transpose() #target_rowのdfを転置
        target_row_dict = target_row_tr.to_dict() #dictのdictとして得られる。keyはindex
        index_number = list(target_row_dict.keys())[0] #indexを数値として取得
        target_row_dict = target_row_dict[index_number] #dictの必要な部分のみ取り出して再代入
        
        window["val"].update(target_row_dict["Object Val"])

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
            window["output_table"].update(composition_df.to_numpy().tolist())
        except:
            pass

    
    if event == "グラフ描画":
            
        num_rows = 4
        num_cols = 5
        
        # ループで10個のサブプロットを生成して表示
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))
        
        for ax, elem in zip(axes.flat,functional_df.columns[2:]):
            x = functional_df[elem].tolist()
            y = functional_df["Object Val"].tolist()
            ax.scatter(x, y, s=3, marker='o', color='blue')  # 仮のデータをプロット
            ax.set_title(elem,fontsize=10)
            ax.set_ylabel('Object Value',fontsize=10)
        # サブプロット間のスペースを調整
        plt.tight_layout()
        plt.get_current_fig_manager().window.wm_geometry("+0+0")
        plt.show(block=False)


    if event == "保存":
        functional_df.to_csv("CSV/functional_groups.csv",index=False)
        composition_df.to_csv("CSV/composition.csv",index=False)


window.close()