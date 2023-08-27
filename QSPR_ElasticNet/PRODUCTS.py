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
sg.theme("DarkTeal10")

layout = [[sg.Text('Products', font=('Constantia',20,"bold"))],
          [sg.Text("")],
          [sg.Button("実行"),
           sg.Button("削除"),
           sg.Button("グラフ描画"),
           sg.Button("Elastic Net")],
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
                  ] for n in range(10)],
                     
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
        
    if event == "実行":

        #出力用dict
        output={} 
        for property in properties:
            calc_val = []
            feed_sum = []
            
            for i in range(10):
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

        functional_group = [elem.replace("(mmol/g)","") for elem in output.keys()][2:]
        quantity = [elem[0] for elem in output.values()][2:]
        
        plt.close()
        fig = plt.figure(figsize=(8, 6)) 
        plt.get_current_fig_manager().window.wm_geometry("+600+0")
        plt.barh(functional_group,quantity,color='gray')
        plt.show(block=False)
        plt.title(output["Name"][0])


        
        #データフレームに新しいデータ行を追加
        new_row = pd.DataFrame(output)
        products_df = pd.concat([products_df,new_row])
        window["output_table"].update([[elem] for elem in products_df["Name"].tolist()])

    
    if event == "削除":
        try:
            products_df = products_df.drop(products_df.index[-1])
            window["output_table"].update([[elem] for elem in products_df["Name"].tolist()])
        except:
            pass

    
    if event == "グラフ描画":
            
        num_rows = 4
        num_cols = 5
        
        # ループで10個のサブプロットを生成して表示
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(12, 8))
        
        for ax, elem in zip(axes.flat,products_df.columns[2:]):
            x = products_df[elem].tolist()
            y = products_df["Object Val"].tolist()
            ax.scatter(x, y, s=3, marker='o', color='blue')  # 仮のデータをプロット
            ax.set_title(elem,fontsize=10)
            ax.set_ylabel('Object Value',fontsize=10)
        # サブプロット間のスペースを調整
        plt.tight_layout()
        plt.get_current_fig_manager().window.wm_geometry("+650+0")
        plt.show(block=False)

    
    if event == "Elastic Net":
        
        y = products_df["Object Val"]
        X = products_df.drop(columns=["Object Val","Name","logP"])
        properties = [word.replace(" (mmol/g)","") for word in X.columns.to_list()]
        
        y_ = y.to_numpy().astype(float)
        X_ = X.to_numpy().astype(float)
        
        from modules import bayes_EN as ben
        ben.Elastic_Net(properties,X_,y_,0.8,100)

products_df.to_csv("CSV/products.csv",index=False)
window.close()