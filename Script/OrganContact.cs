using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class OrganContact : MonoBehaviour
{
    // 静态变量保存类的实例
    public static OrganContact Instance { get; private set; }

    public int num = 0;

    public List<int[]> contact = new List<int[]>(); //0：肝脏点   1：胆囊点

    public int[] contact_liver = {273,340,274,342,277,299,  267,322,318,276,317,326,323,341,324,297,278,298,316,275,
                309,325,268,279,284,344,280,301,327,300,314,285,319,281,265,328,331,303,330,302,329,
                343,266,286,308,282,312,283,339,306,293,307,292,349,315,332,262,290,313,335,
                333,345,287,304,336,347,310,334,320,338,264,337,288,346,289,348,291,305,263,269,
                261,272,311,321,294,270,295,271,296};

    public int[] contact_gallbladder = {392,123,124,236,309,25,126,251,113,150,90,
                                 337,132,137,203,58,229,336,198,76,205,232,
                                 40,93,224,7,283,284,384,252,2,77,169,353,
                                 8,200,74,151,276,263,216,199,152,186,173,
                                 211,212,75,398,381,210,243,35,332,166,373,
                                 361,87,223,259,119,352,282,191,375,371,296,
                                 148,351,313,141,177,295,266,388,395,334,322,
                                 356,189,176,37,394,125,1,272,78,318,248};

    public Vector3[] contact_Pos;

    // Awake在对象初始化时调用
    private void Awake()
    {
        // 检查是否已有实例存在
        if (Instance == null)
        {
            // 如果不存在，则将此对象设置为那个唯一的实例
            Instance = this;
            // 当切换场景时不销毁此对象
            DontDestroyOnLoad(gameObject);
        }
        else
        {
            // 如果已经存在一个实例且不是此实例，销毁此对象
            if (Instance != this)
            {
                Destroy(gameObject);
            }
        }
    }
    private void Start()
    {
        initialize();
        Debug.Log(contact.Count);
    }
    public void initialize()
    {
        num = contact_liver.Length;
        contact_Pos = new Vector3[num];
        for (int i = 0; i < num; ++i)
        {
            int[] x = { contact_liver[i], contact_gallbladder[i] };
            contact.Add(x);
        }
            /*
            liver = new MSM();
            gallbladder = new Low2HighXPBD();

            liver.initialize();
            liverPos = liver.originXpos;
            Debug.Log(liverPos.Length);
            Debug.Log(liverPos[0].x);
            gallbladder.initialize();
            gallbladderPos = gallbladder.High_Xpos;
            Debug.Log(gallbladderPos.Length);
            Debug.Log(gallbladderPos[0].x);

            Debug.Log((liverPos[299] - gallbladderPos[392]).magnitude);


            contact_gallbladder = new int[contact_liver.Length];

            for (int i = 0; i < contact_liver.Length; ++i) {
                Vector3 pos = liverPos[contact_liver[i]];
                for (int j = 0; j < gallbladderPos.Length; ++j) {
                    if (Mathf.Abs((pos - gallbladderPos[j]).magnitude) <= 0.05f) {
                        contact_gallbladder[i] = j;
                        //Debug.Log("--------");
                        continue;
                    }
                }
            }
            hasgall = true;
            */
        }


}
