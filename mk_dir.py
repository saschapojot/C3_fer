from pathlib import Path
from decimal import Decimal, getcontext
from math import factorial
#This script creates directories and conf files for mc

def format_using_decimal(value, precision=25):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

N=5 #unit cell number
TVals=[0.1,1,1.5]

dataRoot="./dataAll/"
dataOutDir=dataRoot

#alpha1: 1/(1!2!)
alpha1_coef=1/(factorial(1)*factorial(2))
alpha1=alpha1_coef
#alpha2: 1/(1!2!)
alpha2_coef=1/(factorial(1)*factorial(2))
alpha2=alpha2_coef
#alpha3: 1/(2!4!)
alpha3_coef=1/(factorial(2)*factorial(4))
alpha3=alpha3_coef
#alpha4: 1/(1!1!2!2!)
alpha4_coef=1/(factorial(1)*factorial(1)*factorial(2)*factorial(2))
alpha4=alpha4_coef
#alpha5: 1/2!

alpha5_coef=1/factorial(2)
alpha5=alpha5_coef
#alpha6: 1/4!
alpha6_coef=1/factorial(4)
alpha6=alpha6_coef
#alpha7: 1/6!
alpha7_coef=1/factorial(6)
alpha7=alpha7_coef

J=1/5*alpha5

print(f"alpha1={alpha1}, alpha2={alpha2}, alpha3={alpha3}, alpha4={alpha4}, alpha5={alpha5}, alpha6={alpha6}, alpha7={alpha7}")
print(f"J={J}")

NStr=format_using_decimal(N)

alpha1_Str=format_using_decimal(alpha1)
alpha2_Str=format_using_decimal(alpha2)
alpha3_Str=format_using_decimal(alpha3)
alpha4_Str=format_using_decimal(alpha4)
alpha5_Str=format_using_decimal(alpha5)
alpha6_Str=format_using_decimal(alpha6)
alpha7_Str=format_using_decimal(alpha7)

J_Str=format_using_decimal(J)

TDirsAll=[]
TStrAll=[]

for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=format_using_decimal(T)
    TStrAll.append(TStr)


def contents_to_conf(k):
    """

    :param k: index of T

    :return:
    """
    contents=[
        "#This is the configuration file for mc computations\n",
        "#System has C3\n",
        "\n" ,
        "#parameters of coefficients\n",
        "\n",
        "#Temperature\n",
        "T="+TStrAll[k]+"\n",
        "\n",
        f"alpha1={alpha1_Str}\n",
        "\n",
        f"alpha2={alpha2_Str}\n",
        "\n",
        f"alpha3={alpha3_Str}\n",
        "\n",
        f"alpha4={alpha4_Str}\n",
        "\n",
        f"alpha5={alpha5_Str}\n",
        "\n",
        f"alpha6={alpha6_Str}\n",
        "\n",
        f"alpha7={alpha7_Str}\n",
        "\n",
        f"J={J_Str}\n",
        "\n",
        f"N={NStr}\n",
        "\n",
        "erase_data_if_exist=False\n",
        "\n",
        "search_and_read_summary_file=True\n"
        "\n",
        "#For the observable name, only digits 0-9, letters a-zA-Z, underscore _ are allowed\n",
        "\n",
        "observable_name=U_dipole\n",
        "\n",
        "effective_data_num_required=1000\n",
        "\n",
        "#this is the data number in each pkl file, i.e., in each flush\n"
        "sweep_to_write=200\n",
        "\n",
        "#within each flush,  sweep_to_write*sweep_multiple mc computations are executed\n",
        "\n",
        "default_flush_num=10\n",
        "\n",
        "h=0.02\n",
        "\n",
        "#the configurations of the system are saved to file if the sweep number is a multiple of sweep_multiple\n",
        "\n",
        "sweep_multiple=300\n",
        ]
    outDir=dataOutDir+f"/N{NStr}/"+"/T"+TStrAll[k]+"/"
    Path(outDir).mkdir(exist_ok=True,parents=True)
    outConfName=outDir+"/run_T"+TStrAll[k]+".mc.conf"
    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)


for k in range(0,len(TVals)):
    contents_to_conf(k)