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
    composition_df = pd.read_csv("CSV/composition.csv")
    

except:
    products_df = pd.DataFrame(columns=["Name","Object Val"]+properties)
    composition_df = pd.DataFrame(columns=["Name","Object Val"])

    

# GUIのレイアウト
sg.theme("DarkTeal10")

layout = [[sg.Text('PRODUCTS', font=('Constantia',20))],
          [sg.Text("")],
          [sg.Button("実行",size = (15,1)),
           sg.Button("データの読み込み",size = (15,1)),
           sg.Button("削除",size = (15,1))],
          [sg.Button("グラフ描画",size = (15,1)),
           sg.Button("Elastic Net",size = (15,1))],
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
        composition ={}

        #組成記録用dictにSample名を記録
        composition["Name"]= values["sample_name"]
        composition["Object Val"]= values["val"]
        
        
        for property in properties:
            calc_val = []
            feed_sum = []
            
            for i in range(10):
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

        functional_group = [elem.replace("(mmol/g)","") for elem in output.keys()][2:]
        quantity = [elem for elem in output.values()][2:]
        
        plt.close()
        fig = plt.figure(figsize=(8, 6)) 
        plt.get_current_fig_manager().window.wm_geometry("+600+0")
        plt.barh(functional_group,quantity,color='gray')
        plt.show(block=False)
        plt.title(output["Name"])

       
        new_name = values["sample_name"]
        #既存データの場合
        if (new_name in products_df["Name"].tolist()) and (new_name in composition_df["Name"].tolist()):
            #読み取ったデータを該当の行に更新
                #官能基
            for key in output.keys():
                products_df.loc[products_df["Name"] == new_name,key] = output[key]
                
                #組成
            for key in composition.keys():
                composition_df.loc[composition_df["Name"] ==new_name,key] = composition[key]

        #新規データの場合
        else: 
            #データフレームに新しいデータ行を追加
            new_row = pd.DataFrame([output])
            new_compositon_row = pd.DataFrame([composition])
            
            #Newrowをproducts_dfに結合、インデックスは振りなおす
            products_df = pd.concat([products_df,new_row],ignore_index=True)            
            composition_df=pd.concat([composition_df,new_compositon_row])

        #各dfの欠損値の0フィリング
        products_df = products_df.fillna(0)
        composition_df = composition_df.fillna(0)

        #記録リストに記録
        window["output_table"].update([[elem] for elem in products_df["Name"].tolist()])


    if event == "データの読み込み":

        edit_sample = values["sample_name"] #検索値をedit_sampleとして格納
        target_row = composition_df.loc[composition_df["Name"] == edit_sample] #検索値と一致する行を取得
        target_row = target_row.loc[:, (target_row != 0).any(axis=0)] #データが０のカラムは削除し、再代入
        target_row_tr = target_row.transpose() #target_rowのdfを転置
        target_row_dict = target_row_tr.to_dict() #dictのdictとして得られる。keyはindex
        index_number = list(target_row_dict.keys())[0] #indexを数値として取得
        target_row_dict = target_row_dict[index_number] #dictの必要な部分のみ取り出して再代入
        
        window["val"].update(target_row_dict["Object Val"])

        #GUIのインプットテーブルを初期化
        for i in range(10):
            window[f"material_{i}"].update("")
            window[f"feed_{i}"].update(0)
        
        #GUIのインプットテーブルの値の更新
        for i, key in enumerate(list(target_row_dict.keys())[2:]):
            window[f"material_{i}"].update(key)
            window[f"feed_{i}"].update(target_row_dict[key])

    
    if event == "削除":
        try:
            #末尾のデータ消去
            products_df = products_df.drop(products_df.index[-1])
            composition_df = composition_df.drop(composition_df.index[-1])
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
composition_df.to_csv("CSV/composition.csv",index=False)

window.close()