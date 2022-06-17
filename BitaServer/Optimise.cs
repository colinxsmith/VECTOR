namespace BitaServer;

public class Optimise
{public Optimise(){
    tlen=0;
    tail=0.05;
    maxRisk=-1;
    minRisk=-1;
    value=-1;
    valuel=-1;
    rmax=-1;
    rmin=-1;
    basket=-1;
    trades=-1;
    min_holding=-1;
    min_trade=-1;
}
public int achievedbasket {get;set;}
public int achievedtrades {get;set;}
public double achievedminhold {get;set;}
public double achievedmintrade {get;set;}

    public double? ogamma { set; get; }
    public double minRisk { set; get; }
    public double maxRisk { set; get; }
    public double rmax { set; get; }
    public double rmin { set; get; }
    public double min_holding { set; get; }
    public double min_trade { set; get; }
    public double? digit { set; get; }
    public double? tdigit { set; get; }
    public int basket { get; set; }
    public int trades { get; set; }
    public int tlen { get; set; }
    public int? nfac { get; set; }
    public int? nabs { get; set; }
    public int? mabs { get; set; }
    public int? n { get; set; }
    public int? m { get; set; }
    public double[]? bench { get; set; }
    public double[]? L { get; set; }
    public double[]? U { get; set; }
    public double[]? A { get; set; }
    public double[]? DATA { get; set; }
    public double[]? buy { get; set; }
    public double[]? Abs_A { get; set; }
    public double[]? Abs_L { get; set; }
    public double[]? Abs_U { get; set; }
    public double[]? mask { get; set; }
    public int?longbasket{set;get;}
    public int?shortbasket{set;get;}
    public int?tradesell{set;get;}
    public int?tradebuy{set;get;}
    public double[]? SV { get; set; }
    public double[]? FC { get; set; }
    public double[]? FL { get; set; }
    public int[]? I_A { get; set; }
    public double[]? sell { get; set; }
    public double? delta { set; get; }
    public double? gamma { set; get; }
    public double value { set; get; }
    public double valuel { set; get; }
    public double? kappa { set; get; }
    public double? Gstrength { get; set; }
    public double tail { get; set; }
    public bool? ETLopt { get; set; }
    public double? ETLmin { get; set; }
    public double? ETLmax { get; set; }
    public string[]? names { get; set; }
    public double[]? initial { get; set; }
    public double[]? Q { get; set; }
    public double[]? w { get; set; }
    public double[]? min_lot { get; set; }
    public double[]? size_lot { get; set; }
    public int[]? shake { get; set; }
    public double? ETL { get; set; }
    public double? VAR { get; set; }
    public int? VARindex { get; set; }
    public double[]? breakdown { get; set; }
    public int? back{set;get;}
    public string? message{set;get;}
    public double?risk{set;get;}
    public double?expreturn{set;get;}
    public double[]?alpha{set;get;}
    public double[]?mctr{set;get;}
    public bool?CVARGLprob{get;set;}
    public double[]? TargetReturn { set; get; }
    public double? LOSS { get; set; }
    public bool? LOSSopt { set; get; }
    public double? LOSSmin { get; set; }
    public double? LOSSmax { set; get; }
    public int?round{get;set;}
}
