import subprocess

# パッケージのリストを指定してpip installを実行する
package_list = ["pubchempy",
                "PySimpleGUI",
                "rdkit",
                "scikit-learn",
                "matplotlib",
                "numpy",
                "optuna",
                "logging",
                "py3Dmol",
                "pandas"]

subprocess.run(['pip', 'install'] + package_list + ['--trusted-host',
                                                    'pypi.python.org',
                                                    '--trusted-host',
                                                    'files.pythonhosted.org',
                                                    '--trusted-host',
                                                    'pypi.org']
                                                    )

import pkg_resources

def is_package_installed(package_name):
    try:
        pkg_resources.get_distribution(package_name)
        return True
    except pkg_resources.DistributionNotFound:
        return False

# パッケージ名を指定して確認する
checker = True

for package_name in  package_list:

    if is_package_installed(package_name):
        print(f"{package_name} finished")
        checker = checker*True
    else:
        print(f"{package_name} error")
        checker = checker*False

if checker ==True:
    print("■■        ■■        ■")
    print(" ■   ■■   ■         ■")
    print(" ■   ■■   ■   ■■■   ■    ■■■■  ■■■■   ■■■■■■■■    ■■■")
    print(" ■   ■■   ■  ■■  ■  ■   ■■    ■■  ■   ■   ■   ■  ■■  ■")
    print(" ■■ ■■ ■  ■  ■   ■■ ■   ■     ■    ■  ■   ■   ■  ■   ■■")
    print(" ■■ ■  ■ ■■  ■■■■■■ ■   ■     ■    ■  ■   ■   ■  ■■■■■■")
    print("  ■ ■  ■ ■   ■      ■   ■     ■    ■  ■   ■   ■  ■     ")
    print("  ■■    ■■   ■■     ■   ■■    ■■  ■   ■   ■   ■  ■■    ")
    print("  ■■    ■■    ■■■■   ■   ■■■■  ■■■■   ■   ■   ■   ■■■■ ")
import time
# 5秒間待機
time.sleep(10)