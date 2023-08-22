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
    products_df = pd.read_csv("CSV/products.csv")
except:
    products_df = pd.DataFrame(columns=["Name","Object Val"]+properties)

# GUIのレイアウト
sg.theme("DarkBrown")

layout = [[sg.Text('PRODUCTS', font=('Constantia',20))],
          [sg.Text("")],
          [sg.Button("Calculate"),sg.Button("Delete")],
          [sg.Text("Sample Name",size=(20, 1)), sg.Input(size=(20, 1),key = "sample_name" )],
          [sg.Text("Value",size=(20, 1)), sg.Input(size=(20, 1),key = "val" )],
          [sg.Text('_'  * 70)], #横線区切り
          [sg.Column(layout=[
                  [sg.Combo(materials, size=(20, 1), key=f'material_{n}'),
                   sg.Input("0", size=(20, 1), key=f"feed_{n}")] for n in range(15)],
              size=(400, 400)  # 列全体のサイズ
              ),
           sg.Table(headings =["Name"],values = [[elem] for elem in products_df["Name"].tolist()],
                    key ="output_table",
                    size=(10,25))
          ]]
# ウインドウの出現位置を指定
win_location = (0, 0)

# モニターの解像度を取得
screen_width, screen_height = sg.Window.get_screen_size()

# ウィンドウのサイズをモニターの半分に設定
window_size = (screen_width // 2, screen_height)

window = sg.Window("PRODUCTS", 
                   layout, 
                   size=window_size,
                   resizable=True,
                   location=win_location)



while True:
    
    event, values = window.read()
    
    if event == sg.WINDOW_CLOSED:
        break
        
    if event == "Calculate":

        #出力用dict
        output={} 
        for property in properties:
            calc_val = []
            feed_sum = []
            
            for i in range(15):
                material_name = values[f"material_{i}"]
                if material_name != "":  #materialnameに入力されている場合のみデータ参照可能
                    feed = float(values[f"feed_{i}"])
                    calc_val.append(df.loc[material_name,property]*feed)
                    feed_sum.append(feed)

            #目的変数値の記録
            output["Name"] = [values["sample_name"]]
            #目的変数値の記録
            output["Object Val"] = [values["val"]]
            #官能基情報の記録
            output[property] = [sum(calc_val)/sum(feed_sum)]

        x = [elem.replace("(mmol/g)","") for elem in output.keys()][2:]
        y = [elem[0] for elem in output.values()][2:]
        
        plt.close()
        fig = plt.figure(figsize=(8, 6)) 
        plt.get_current_fig_manager().window.wm_geometry("+550+0")
        plt.barh(x,y,color='gray')
        plt.show(block=False)
        plt.title(output["Name"][0])


        
        #データフレームに新しいデータ行を追加
        new_row = pd.DataFrame(output)
        products_df = pd.concat([products_df,new_row])
        window["output_table"].update([[elem] for elem in products_df["Name"].tolist()])

    
    if event == "Delete":
        try:
            products_df = products_df.drop(products_df.index[-1])
            window["output_table"].update([[elem] for elem in products_df["Name"].tolist()])
        except:
            pass

products_df.to_csv("CSV/products.csv",index=False)
window.close()