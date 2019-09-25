unit UnitBlockViewFactors;
// К математико-вычислительной модели излучения внутри vacuum Prism Object.
// Визуализирует и рассчитывает View Factors для Prism Object
// на основе априорно аналитически взятых значений интеграла.
// Heat and Mass Transfer: Fundamentals & Applications. Fourth Edition.
// Yunus A. Cengel, Afshin J.Ghajar. McGraw-Hill, 2011.
// Isidoro Martinez Radiative View Factors. 1995-2016.
// Также в данном модуле находятся используя заданные :
// 1. emissivity на гранях Prism.
// 2. Температуры в К на гранях Prism.
// 3. Вычисленные ранее View Factors на гранях Prism.
// находтся плотности тепловых потоков JW, JE, JS, JN, JB, JT излучаемых с каждой грани
// внутрь вакуумного Prism элемента.
// Во время сборки матрицы для уравнения теплопередачи тепловые потоки JG, G: W,E,S,N,B,T
// аддитивно добавляются  с учётом направления к кондукционных и конвективным потокам
// на гранях Prism элемента. Внутри AliceFlow_v0_24.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ExtCtrls, Vcl.StdCtrls;

type
  TFormRadiation = class(TForm)
    GroupBox1: TGroupBox;
    RadioGroup1: TRadioGroup;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    LabelSum: TLabel;
    GroupBoxemissivity: TGroupBox;
    Label6: TLabel;
    Label7: TLabel;
    Edit1: TEdit;
    Edit2: TEdit;
    Label8: TLabel;
    Label9: TLabel;
    Edit3: TEdit;
    Edit4: TEdit;
    Label10: TLabel;
    Label11: TLabel;
    Edit5: TEdit;
    Edit6: TEdit;
    GroupBox2: TGroupBox;
    Label12: TLabel;
    Label13: TLabel;
    Edit7: TEdit;
    Label14: TLabel;
    Edit8: TEdit;
    Label15: TLabel;
    Label16: TLabel;
    Label17: TLabel;
    Edit9: TEdit;
    Edit10: TEdit;
    Label18: TLabel;
    Label19: TLabel;
    Edit11: TEdit;
    Edit12: TEdit;
    GroupBox3: TGroupBox;
    Label20: TLabel;
    Label21: TLabel;
    Label22: TLabel;
    Label23: TLabel;
    Label24: TLabel;
    Label25: TLabel;
    Label26: TLabel;
    Label27: TLabel;
    Label28: TLabel;
    Label29: TLabel;
    Button1: TButton;
    Editrelaxfactor: TEdit;
    Label30: TLabel;
    CheckBoxinit: TCheckBox;
    ButtonApply: TButton;
    CheckBoxinternalRadiation: TCheckBox;
    procedure RadioGroup1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure CheckBoxinternalRadiationClick(Sender: TObject);
    procedure ButtonApplyClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
    k1 : Integer;
     // The density heat flux
    JW, JE, JS, JN, JB, JT : Real;
    // Вычисляет view factor между двумя параллельными пластинами.
    function view_factor_aligned_parallel_rectangles(x34 : Extended; y34 : Extended; L34 : Extended) : Extended;
    // Вычисляет view factor между двумя пластинами составляющими между собой прямой угол
    // и контактирующих по общему ребру.
    function view_factor_perpendicular_rectangles_with_a_common_edge(x34 : Extended; y34 : Extended; z34 : Extended) : Extended;
  end;

var
  FormRadiation: TFormRadiation;

implementation

{$R *.dfm}
uses
    Math, VisualUnit,  Unitrectangularplot, UnitVariables;

// Находит плотности тепловых потоков излучения  на гранях Prism Object при
// 1. Заданных температурах в К на гранях Prism Object.
// 2. Излучательных способностях на гранях Prism Object.
// 3. Вычисленных значениях View Factors.
procedure TFormRadiation.Button1Click(Sender: TObject);
var
    bOk, bOk2, bOk3 : Boolean;
    s : String;
    r, alpha : Real;
    k : Integer;
    x34, y34, z34 : Real;
    // The View Factors.
    FWE, FWS, FWN, FWB, FWT : Real;
    FEW, FES, FEN, FEB, FET : Real;
    FSW, FSE, FSN, FSB, FST : Real;
    FNW, FNE, FNS, FNB, FNT : Real;
    FBW, FBE, FBS, FBN, FBT : Real;
    FTW, FTE, FTS, FTN, FTB : Real;
    // The emissivities
    epsW, epsE, epsS, epsN, epsB, epsT : Real;
    // The Temperatures
    TW, TE, TS, TN, TB, TT : Real;
    // The Stefan-Bolcman const.
    sigma : Real;
    // The density heat flux
    JW1, JE1, JS1, JN1, JB1, JT1 : Real;
    // Диагональные коэффициенты матрицы.
    apW, apE, apS, apN, apB, apT : Real;


begin
   // Проверка корректности ввода.
   s:=Trim(Edit1.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit1.Text:=s;

   s:=Trim(Edit2.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit2.Text:=s;
   s:=Trim(Edit3.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit3.Text:=s;
   s:=Trim(Edit4.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit4.Text:=s;
   s:=Trim(Edit5.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit5.Text:=s;
   s:=Trim(Edit6.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit6.Text:=s;
   s:=Trim(Edit7.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit7.Text:=s;
   s:=Trim(Edit8.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit8.Text:=s;
   s:=Trim(Edit9.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit9.Text:=s;
   s:=Trim(Edit10.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit10.Text:=s;
   s:=Trim(Edit11.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit11.Text:=s;
   s:=Trim(Editrelaxfactor.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Editrelaxfactor.Text:=s;


   // Вычисление all View Factors:
   // min X (W)
   x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS); // x
   y34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS); // y
   z34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS); // z
   FWE:=view_factor_aligned_parallel_rectangles(x34,y34,z34);
   FWS:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FWN:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FWB:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FWT:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   // max X (E)
   x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
   y34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
   z34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
   FEW:=view_factor_aligned_parallel_rectangles(y34,x34,z34);
   FES:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FEN:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FEB:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FET:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   // min Y (S)
   z34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
   x34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
   y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
   FSW:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FSE:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FSN:=view_factor_aligned_parallel_rectangles(x34,y34,z34);
   FSB:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FST:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   // max Y (N)
   z34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
   x34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
   y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
   FNW:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FNE:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FNS:=view_factor_aligned_parallel_rectangles(x34,y34,z34);
   FNB:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FNT:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   // min Z (B)
   x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
   z34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
   y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
   FBW:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FBE:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FBS:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FBN:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FBT:=view_factor_aligned_parallel_rectangles(x34,y34,z34);
   // max Z (T)
   x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
   z34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
   y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
   FTW:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FTE:=view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34);
   FTS:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FTN:=view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34);
   FTB:=view_factor_aligned_parallel_rectangles(x34,y34,z34);

   // Считывание значений emissivity:
    //epsW:=StrToFloat(Edit1.Text);
    //epsE:=StrToFloat(Edit2.Text);
    //epsS:=StrToFloat(Edit3.Text);
    //epsN:=StrToFloat(Edit4.Text);
    //epsB:=StrToFloat(Edit5.Text);
    //epsT:=StrToFloat(Edit6.Text);
    bOk3:=true;
    if bOk3 then epsW:=FormVariables.my_real_convert(Edit1.Text,bOk3);
    if bOk3 then epsE:=FormVariables.my_real_convert(Edit2.Text,bOk3);
    if bOk3 then epsS:=FormVariables.my_real_convert(Edit3.Text,bOk3);
    if bOk3 then epsN:=FormVariables.my_real_convert(Edit4.Text,bOk3);
    if bOk3 then epsB:=FormVariables.my_real_convert(Edit5.Text,bOk3);
    if bOk3 then epsT:=FormVariables.my_real_convert(Edit6.Text,bOk3);

    //r:=StrToFloat(Edit1.Text);
    r:=FormVariables.my_real_convert(Edit1.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit1.Color:=clred;
    end
     else
    begin
       Edit1.Color:=clWhite;
    end;
     //r:=StrToFloat(Edit2.Text);
    r:=FormVariables.my_real_convert(Edit2.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit2.Color:=clred;
    end
     else
    begin
       Edit2.Color:=clWhite;
    end;
     //r:=StrToFloat(Edit3.Text);
     r:=FormVariables.my_real_convert(Edit3.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit3.Color:=clred;
    end
     else
    begin
       Edit3.Color:=clWhite;
    end;
    //r:=StrToFloat(Edit4.Text);
    r:=FormVariables.my_real_convert(Edit4.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit4.Color:=clred;
    end
     else
    begin
       Edit4.Color:=clWhite;
    end;
     //r:=StrToFloat(Edit5.Text);
     r:=FormVariables.my_real_convert(Edit5.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit5.Color:=clred;
    end
     else
    begin
       Edit5.Color:=clWhite;
    end;
     //r:=StrToFloat(Edit6.Text);
     r:=FormVariables.my_real_convert(Edit6.Text,bOk3);
    if (not(bOk3)or(r>1.0)or(r<0.0)) then
    begin
       bOk2:=False;
       Edit6.Color:=clred;
    end
     else
    begin
       Edit6.Color:=clWhite;
    end;

    sigma:=5.670367e-8; // Wxm!(-2)xK!(-4).
    bOk2:=true;
    // Считывание температур :
    TW:=StrToFloat(Edit7.Text);
    TE:=StrToFloat(Edit8.Text);
    TS:=StrToFloat(Edit9.Text);
    TN:=StrToFloat(Edit10.Text);
    TB:=StrToFloat(Edit11.Text);
    TT:=StrToFloat(Edit12.Text);

    if (TW<0.0) then
    begin
       bOk2:=false;
       Edit7.Color:=clred;
    end
    else
    begin
      Edit7.Color:=clWhite;
    end;
    if (TE<0.0) then
    begin
        bOk2:=false;
        Edit8.Color:=clred;
    end
    else
    begin
         Edit8.Color:=clWhite;
    end;
    if (TS<0.0) then
    begin
       bOk2:=false;
       Edit9.Color:=clred;
    end
     else
    begin
        Edit9.Color:=clWhite;
    end;
    if (TN<0.0) then
    begin
       bOk2:=false;
       Edit10.Color:=clred;
    end
     else
    begin
        Edit10.Color:=clWhite;
    end;
    if (TB<0.0) then
    begin
       bOk2:=false;
       Edit11.Color:=clred;
    end
     else
    begin
        Edit11.Color:=clWhite;
    end;
    if (TT<0.0) then
    begin
       bOk2:=false;
       Edit12.Color:=clred;
    end
     else
    begin
        Edit12.Color:=clWhite;
    end;

    if (bOk2=true) then
    begin
       // The density heat flux
       // initialization:
       if (CheckBoxInit.Checked=true) then
       begin
          JW:=0.0;
          JE:=0.0;
          JS:=0.0;
          JN:=0.0;
          JB:=0.0;
          JT:=0.0;


          frmRectangularPlot.Show;
          frmRectangularPlot.cht1.Title.Text[0]:='Radiation density heat flux residual';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='Evklid residual';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='iteration';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Axes.Left.LogarithmicBase:=10.0;
          frmRectangularPlot.cht1.Axes.Left.Logarithmic:=True;

          k1:=0;
          bOk:=true;
       end;

       alpha:=StrToFloat(Editrelaxfactor.Text);
       r:=StrToFloat(Editrelaxfactor.Text);
        if ((r>1.0)or(r<0.0)) then
        begin
           bOk2:=False;
           Editrelaxfactor.Color:=clred;
           alpha:=0.1;
           ShowMessage('set alpha:=0.1 default value. alpha must be 0.0<alpha<1.0');
        end
         else
        begin
           Editrelaxfactor.Color:=clWhite;
        end;
       for k := 0 to 100 do
       begin

          //JW1:=alpha*sigma*TW*TW*TW*TW-((1.0-epsW)/epsW)*(FWE*(JW-JE)+FWS*(JW-JS)+FWN*(JW-JN)+FWB*(JW-JB)+FWT*(JW-JT))+(1.0-alpha)*JW;
          //  JE1:=alpha*sigma*TE*TE*TE*TE-((1.0-epsE)/epsE)*(FEW*(JE-JW)+FES*(JE-JS)+FEN*(JE-JN)+FEB*(JE-JB)+FET*(JE-JT))+(1.0-alpha)*JE;
          //JS1:=alpha*sigma*TS*TS*TS*TS-((1.0-epsS)/epsS)*(FSW*(JS-JW)+FSE*(JS-JE)+FSN*(JS-JN)+FSB*(JS-JB)+FST*(JS-JT))+(1.0-alpha)*JS;
          //JN1:=alpha*sigma*TN*TN*TN*TN-((1.0-epsN)/epsN)*(FNW*(JN-JW)+FNE*(JN-JE)+FNS*(JN-JS)+FNB*(JN-JB)+FNT*(JN-JT))+(1.0-alpha)*JN;
          //JB1:=alpha*sigma*TB*TB*TB*TB-((1.0-epsB)/epsB)*(FBW*(JB-JW)+FBE*(JB-JE)+FBS*(JB-JS)+FBN*(JB-JN)+FBT*(JB-JT))+(1.0-alpha)*JB;
          //JT1:=alpha*sigma*TT*TT*TT*TT-((1.0-epsT)/epsT)*(FTW*(JT-JW)+FTE*(JT-JE)+FTS*(JT-JS)+FTN*(JT-JN)+FTB*(JT-JB))+(1.0-alpha)*JT;

          // Вычисляем диагональные коэффициенты.
          // apW:=1.0+((1.0-epsE)/epsE)*(FWE+FWS+FWN+FWB+FWT);
          // simplify
          apW:=1.0+ ((1.0-epsW)/epsW)*1.0;
          apE:=1.0+ ((1.0-epsE)/epsE)*1.0;
          apS:=1.0+ ((1.0-epsS)/epsS)*1.0;
          apN:=1.0+ ((1.0-epsN)/epsN)*1.0;
          apB:=1.0+ ((1.0-epsB)/epsB)*1.0;
          apT:=1.0+ ((1.0-epsT)/epsT)*1.0;


          JW1:=alpha*((sigma*TW*TW*TW*TW+((1.0-epsW)/epsW)*(FWE*(JE)+FWS*(JS)+FWN*(JN)+FWB*(JB)+FWT*(JT)))/apW)+(1.0-alpha)*JW;
          JE1:=alpha*((sigma*TE*TE*TE*TE+((1.0-epsE)/epsE)*(FEW*(JW)+FES*(JS)+FEN*(JN)+FEB*(JB)+FET*(JT)))/apE)+(1.0-alpha)*JE;
          JS1:=alpha*((sigma*TS*TS*TS*TS+((1.0-epsS)/epsS)*(FSW*(JW)+FSE*(JE)+FSN*(JN)+FSB*(JB)+FST*(JT)))/apS)+(1.0-alpha)*JS;
          JN1:=alpha*((sigma*TN*TN*TN*TN+((1.0-epsN)/epsN)*(FNW*(JW)+FNE*(JE)+FNS*(JS)+FNB*(JB)+FNT*(JT)))/apN)+(1.0-alpha)*JN;
          JB1:=alpha*((sigma*TB*TB*TB*TB+((1.0-epsB)/epsB)*(FBW*(JW)+FBE*(JE)+FBS*(JS)+FBN*(JN)+FBT*(JT)))/apB)+(1.0-alpha)*JB;
          JT1:=alpha*((sigma*TT*TT*TT*TT+((1.0-epsT)/epsT)*(FTW*(JW)+FTE*(JE)+FTS*(JS)+FTN*(JN)+FTB*(JB)))/apT)+(1.0-alpha)*JT;



          // невязка :
          r:=Sqrt(Sqr(JW1-sigma*TW*TW*TW*TW-((1.0-epsW)/epsW)*(FWE*(JW-JE)+FWS*(JW-JS)+FWN*(JW-JN)+FWB*(JW-JB)+FWT*(JW-JT)))+
          Sqr(JE1-sigma*TE*TE*TE*TE-((1.0-epsE)/epsE)*(FEW*(JE-JW)+FES*(JE-JS)+FEN*(JE-JN)+FEB*(JE-JB)+FET*(JE-JT)))+
          Sqr(JS1-sigma*TS*TS*TS*TS-((1.0-epsS)/epsS)*(FSW*(JS-JW)+FSE*(JS-JE)+FSN*(JS-JN)+FSB*(JS-JB)+FST*(JS-JT)))+
          Sqr(JN1-sigma*TN*TN*TN*TN-((1.0-epsN)/epsN)*(FNW*(JN-JW)+FNE*(JN-JE)+FNS*(JN-JS)+FNB*(JN-JB)+FNT*(JN-JT)))+
          Sqr(JB1-sigma*TB*TB*TB*TB-((1.0-epsB)/epsB)*(FBW*(JB-JW)+FBE*(JB-JE)+FBS*(JB-JS)+FBN*(JB-JN)+FBT*(JB-JT)))+
          Sqr(JT1-sigma*TT*TT*TT*TT-((1.0-epsT)/epsT)*(FTW*(JT-JW)+FTE*(JT-JE)+FTS*(JT-JS)+FTN*(JT-JN)+FTB*(JT-JB))));
          frmRectangularPlot.cht1.Series[0].AddXY(k+k1,r,'',clRed);
          if (bOk) then
          begin
             frmRectangularPlot.cht1.Axes.Left.Maximum:=1.0e6;
             bOk:=false;
          end;
          frmRectangularPlot.cht1.Axes.Left.Minimum:=1.0e-10;

          JW:=JW1;
          JE:=JE1;
          JS:=JS1;
          JN:=JN1;
          JB:=JB1;
          JT:=JT1;
        end;
        k1:=k1+101;

        label24.Caption:='JW='+FloatToStr(JW);
        label25.Caption:='JE='+FloatToStr(JE);
        label26.Caption:='JS='+FloatToStr(JS);
        label27.Caption:='JN='+FloatToStr(JN);
        label28.Caption:='JB='+FloatToStr(JB);
        label29.Caption:='JT='+FloatToStr(JT);
    end
    else
    begin
       ShowMessage('Your input data is incorrect. Please, reinput...');
    end;

end;

// Ввод значений emissivity.
procedure TFormRadiation.ButtonApplyClick(Sender: TObject);
var
   s : String;
   k : Integer;
   bOk : Boolean;
   r : Real;
begin
   // Ввод значений emissivity.
    // Проверка корректности ввода.
   s:=Trim(Edit1.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit1.Text:=s;

   s:=Trim(Edit2.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit2.Text:=s;
   s:=Trim(Edit3.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit3.Text:=s;
   s:=Trim(Edit4.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit4.Text:=s;
   s:=Trim(Edit5.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit5.Text:=s;
   s:=Trim(Edit6.Text);
   for k:=1 to length(s) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s[k]='.') then s[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s[k]=',') then s[k]:='.';
      end;
   end;
   Edit6.Text:=s;

   bOk:=true;

   //r:=StrToFloat(Edit1.Text);
   r:=FormVariables.my_real_convert(Edit1.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit1.Color:=clred;
   end
   else
   begin
      Edit1.Color:=clWhite;
   end;

    //r:=StrToFloat(Edit2.Text);
    r:=FormVariables.my_real_convert(Edit2.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit2.Color:=clred;
   end
   else
   begin
      Edit2.Color:=clWhite;
   end;

   //r:=StrToFloat(Edit3.Text);
   r:=FormVariables.my_real_convert(Edit3.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit3.Color:=clred;
   end
   else
   begin
      Edit3.Color:=clWhite;
   end;

    //r:=StrToFloat(Edit4.Text);
    r:=FormVariables.my_real_convert(Edit4.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit4.Color:=clred;
   end
   else
   begin
      Edit4.Color:=clWhite;
   end;

   //r:=StrToFloat(Edit5.Text);
   r:=FormVariables.my_real_convert(Edit5.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit5.Color:=clred;
   end
   else
   begin
      Edit5.Color:=clWhite;
   end;

   //r:=StrToFloat(Edit6.Text);
   r:=FormVariables.my_real_convert(Edit6.Text,bOk);
   if ((r>1.0)or(r<0.0)) then
   begin
     bOk:=False;
     Edit6.Color:=clred;
   end
   else
   begin
      Edit6.Color:=clWhite;
   end;

   if (bOk) then
   begin
       // Считывание значений emissivity:
       Laplas.body[Laplas.itek].emissW:=FormVariables.my_real_convert(Edit1.Text,bOk);
       Laplas.body[Laplas.itek].emissE:=FormVariables.my_real_convert(Edit2.Text,bOk);
       Laplas.body[Laplas.itek].emissS:=FormVariables.my_real_convert(Edit3.Text,bOk);
       Laplas.body[Laplas.itek].emissN:=FormVariables.my_real_convert(Edit4.Text,bOk);
       Laplas.body[Laplas.itek].emissB:=FormVariables.my_real_convert(Edit5.Text,bOk);
       Laplas.body[Laplas.itek].emissT:=FormVariables.my_real_convert(Edit6.Text,bOk);

       // излучательная способность
       Laplas.body[Laplas.itek].semissW:=Trim(Edit1.Text);
       Laplas.body[Laplas.itek].semissE:=Trim(Edit2.Text);
       Laplas.body[Laplas.itek].semissS:=Trim(Edit3.Text);
       Laplas.body[Laplas.itek].semissN:=Trim(Edit4.Text);
       Laplas.body[Laplas.itek].semissB:=Trim(Edit5.Text);
       Laplas.body[Laplas.itek].semissT:=Trim(Edit6.Text);

       if (abs(Laplas.body[Laplas.itek].emissW)<1.0e-12) then
       begin
          // Во избежании деления на ноль.
          Laplas.body[Laplas.itek].emissW:=0.001;
            if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissW:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissW:='0.001';
         end;
       end;
       if (abs(Laplas.body[Laplas.itek].emissE)<1.0e-12) then
       begin
          // Во избежании деления на ноль.
          Laplas.body[Laplas.itek].emissE:=0.001;
          if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissE:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissE:='0.001';
         end;
       end;
       if (abs(Laplas.body[Laplas.itek].emissS)<1.0e-12) then
       begin
         // Во избежании деления на ноль.
         Laplas.body[Laplas.itek].emissS:=0.001;
         if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissS:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissS:='0.001';
         end;
       end;
       if (abs(Laplas.body[Laplas.itek].emissN)<1.0e-12) then
       begin
          // Во избежании деления на ноль.
          Laplas.body[Laplas.itek].emissN:=0.001;
          if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissN:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissN:='0.001';
         end;
       end;
       if (abs(Laplas.body[Laplas.itek].emissB)<1.0e-12) then
       begin
          // Во избежании деления на ноль.
          Laplas.body[Laplas.itek].emissB:=0.001;
          if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissB:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissB:='0.001';
         end;
       end;
       if (abs(Laplas.body[Laplas.itek].emissT)<1.0e-12) then
       begin
          // Во избежании деления на ноль.
          Laplas.body[Laplas.itek].emissT:=0.001;
          if (FormatSettings.DecimalSeparator=',') then
          begin
             Laplas.body[Laplas.itek].semissT:='0,001';
          end;
          if (FormatSettings.DecimalSeparator='.') then
          begin
             Laplas.body[Laplas.itek].semissT:='0.001';
         end;
       end;

       if (CheckBoxinternalRadiation.Checked=true) then
       begin
          // Учитывается модель радиационного теплообмена
          // внутри Prism Object.
          Laplas.body[Laplas.itek].binternalRadiation:=1;
       end
        else
       begin
         // Не учитывается модель радиационного теплообмена
         // внутри Prism Object.
         Laplas.body[Laplas.itek].binternalRadiation:=0;
       end;
    end
     else
    begin
       ShowMessage('Your input is incorrect. Please, reinput...');
    end;

end;

// Показывает пользователю вычисленные значения ViewFactors.
procedure TFormRadiation.CheckBoxinternalRadiationClick(Sender: TObject);
begin
   if (CheckBoxInternalRadiation.Checked=true) then
   begin
       GroupBox1.Visible:=true;
       GroupBox2.Visible:=true;
       GroupBox3.Visible:=true;
   end
    else
   begin
       GroupBox1.Visible:=false;
       GroupBox2.Visible:=false;
       GroupBox3.Visible:=false;
   end;
end;

procedure TFormRadiation.FormCreate(Sender: TObject);
begin
   // значение emissivity по умолчанию.
   // Отражательная пспособность равна rho:=1.0-emissivity.
   if (FormatSettings.DecimalSeparator=',') then
   begin
      Edit1.Text:='0,8';
      Edit2.Text:='0,8';
      Edit3.Text:='0,8';
      Edit4.Text:='0,8';
      Edit5.Text:='0,8';
      Edit6.Text:='0,8';
   end;
   if (FormatSettings.DecimalSeparator='.') then
   begin
       Edit1.Text:='0.8';
       Edit2.Text:='0.8';
       Edit3.Text:='0.8';
       Edit4.Text:='0.8';
       Edit5.Text:='0.8';
       Edit6.Text:='0.8';
   end;
   k1:=0;
end;



procedure TFormRadiation.RadioGroup1Click(Sender: TObject);
var
   x34,y34,z34,s34 : Extended;
begin
   case RadioGroup1.ItemIndex of
      0 : begin
             // minX (W)
             x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS); // x
             y34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS); // y
             z34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS); // z
             label1.Caption:='FWE='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(x34,y34,z34));
             label2.Caption:='FWS='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label3.Caption:='FWN='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label4.Caption:='FWB='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label5.Caption:='FWT='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             s34:= view_factor_aligned_parallel_rectangles(x34,y34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
      1 : begin
             // MaxX (E)
             x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
             y34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
             z34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
             label1.Caption:='FEW='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(y34,x34,z34));
             label2.Caption:='FES='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label3.Caption:='FEN='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label4.Caption:='FEB='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label5.Caption:='FET='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             s34:=view_factor_aligned_parallel_rectangles(y34,x34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
      2 : begin
              // MinY (S)
             z34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
             x34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
             y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
             label1.Caption:='FSW='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label2.Caption:='FSE='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label3.Caption:='FSN='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(x34,y34,z34));
             label4.Caption:='FSB='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label5.Caption:='FST='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             s34:= view_factor_aligned_parallel_rectangles(x34,y34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
      3 : begin
             // MaxY (N)
             z34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
             x34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
             y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
             label1.Caption:='FNW='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label2.Caption:='FNE='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label3.Caption:='FNS='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(x34,y34,z34));
             label4.Caption:='FNB='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label5.Caption:='FNT='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             s34:=view_factor_aligned_parallel_rectangles(x34,y34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
      4 : begin
             // MinZ (B)
             x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
             z34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
             y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
             label1.Caption:='FBW='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label2.Caption:='FBE='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label3.Caption:='FBS='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label4.Caption:='FBN='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label5.Caption:='FBT='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(x34,y34,z34));
             s34:=view_factor_aligned_parallel_rectangles(x34,y34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
      5 : begin
             // MaxZ (T)
             x34:=Abs(Laplas.body[Laplas.itek].yE-Laplas.body[Laplas.itek].yS);
             z34:=Abs(Laplas.body[Laplas.itek].zE-Laplas.body[Laplas.itek].zS);
             y34:=Abs(Laplas.body[Laplas.itek].xE-Laplas.body[Laplas.itek].xS);
             label1.Caption:='FTW='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label2.Caption:='FTE='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             label3.Caption:='FTS='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label4.Caption:='FTN='+FormatFloat('#.###',view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34));
             label5.Caption:='FTB='+FormatFloat('#.###',view_factor_aligned_parallel_rectangles(x34,y34,z34));
             s34:=view_factor_aligned_parallel_rectangles(x34,y34,z34)+2.0*(view_factor_perpendicular_rectangles_with_a_common_edge(y34,x34,z34)+view_factor_perpendicular_rectangles_with_a_common_edge(x34,y34,z34));
             labelSum.Caption:='validate='+FormatFloat('#.###',s34);
          end;
   end;
end;

// Вычисляет view factor между двумя параллельными пластинами.
function TFormRadiation.view_factor_aligned_parallel_rectangles(x34 : Extended; y34 : Extended; L34 : Extended) : Extended;
var
  x_tilda, y_tilda : Extended;
  x_tilda1, y_tilda1 : Extended;
  m_PI : Extended; // число ПИ.
  Fi_j : Extended; // View Factor
begin
    // Вычисление View Factors между двумя паралельными пластинами на основе аналитически
    // вычисленного значения интеграла.
    // Heat and Mass Transfer: Fundamentals & Applications. Fourth Edition.
    // Yunus A. Cengel, Afshin J.Ghajar. McGraw-Hill, 2011.
    x_tilda:=x34/L34;
    y_tilda:=y34/L34;
    if ((x_tilda>120.0)or (y_tilda>120.0)) then
    begin
       ShowMessage('X='+FloatToStr(x_tilda)+' Y='+FloatToStr(y_tilda));
    end;
    m_PI:=3.1415926;
    // неверная формула. (БРЕД)
   // Fi_j:=(2.0/(m_PI*x_tilda*y_tilda))*
    //((Ln(Sqrt(((1.0+x_tilda*x_tilda)/(1.0+x_tilda*x_tilda+y_tilda*y_tilda))*
    //(1.0+y_tilda*y_tilda))))
    //+(x_tilda*Sqrt(1.0+y_tilda*y_tilda))/(Tan(x_tilda/Sqrt(1.0+y_tilda*y_tilda)))
    //+(y_tilda*Sqrt(1.0+x_tilda*x_tilda))/(Tan(y_tilda/Sqrt(1.0+x_tilda*x_tilda)))
    //-x_tilda/Tan(x_tilda)-y_tilda/Tan(y_tilda));
    // Isidoro Martinez Radiative View Factors. 1995-2016.
    x_tilda1:=sqrt(1.0+x_tilda*x_tilda);
    y_tilda1:=sqrt(1.0+y_tilda*y_tilda);
    Fi_j:=(1.0/(m_PI*x_tilda*y_tilda))*(Ln((x_tilda1*x_tilda1*y_tilda1*y_tilda1)/
    (x_tilda1*x_tilda1+y_tilda1*y_tilda1-1.0))+
    +2.0*x_tilda*(y_tilda1*ArcTan(x_tilda/y_tilda1)-ArcTan(x_tilda))+
    2.0*y_tilda*(x_tilda1*ArcTan(y_tilda/x_tilda1)-ArcTan(y_tilda)));
    Result:=Fi_j;
end;


// Вычисляет view factor между двумя пластинами составляющими между собой прямой угол
// и контактирующих по общему ребру.
function TFormRadiation.view_factor_perpendicular_rectangles_with_a_common_edge(x34 : Extended;
 y34 : Extended; z34 : Extended) : Extended;
 var
    h, w,a,b,c : Extended;
    m_PI : Extended; // число ПИ.
    Fi_j : Extended; // View Factor
 begin
    // Вычисление View Factors между двумя перпендикулярными пластинами имеющими общее ребро
    // на основе аналитически
    // вычисленного значения интеграла.
    // Heat and Mass Transfer: Fundamentals & Applications. Fourth Edition.
    // Yunus A. Cengel, Afshin J.Ghajar. McGraw-Hill, 2011.
    h:=z34/x34;
    w:=y34/x34;
   if ((h>140.0)or (w>140.0)) then
    begin
       ShowMessage('H='+FloatToStr(h)+' W='+FloatToStr(w));
    end;
    m_PI:=3.1415926;
    // неверная формула (БРЕД)
    //Fi_j:=(1.0/(m_PI*w))*(w/Tan(1.0/w)+h/Tan(1.0/h)-Sqrt(h*h+w*w)/
    //Tan(1.0/(Sqrt(h*h+w*w)))
    //+0.25*Ln((((1.0+w*w)*(1+h*h))/(1.0+w*w+h*h))
    //*Power((((w*w)/(1+w*w))*((1.0+w*w+h*h)/(w*w+h*h))),(w*w))
    //*Power((((h*h)/(1.0+h*h))*((1.0+h*h+w*w)/(h*h+w*w))),(h*h))));
    // Isidoro Martinez Radiative View Factors. 1995-2016.
    a:=((1.0+h*h)*(1.0+w*w))/(1.0+h*h+w*w);
    b:=(w*w*(1.0+h*h+w*w))/((1.0+w*w)*(h*h+w*w));
    c:=(h*h*(1.0+h*h+w*w))/((1.0+h*h)*(h*h+w*w));
    Fi_j:=(1.0/(m_PI*w))*(h*ArcTan(1.0/h)+w*ArcTan(1.0/w)-Sqrt(h*h+w*w)*
    ArcTan(1.0/(Sqrt(h*h+w*w)))+
    0.25*Ln(a*Power(b,w*w)*Power(c,h*h)));
    Result:=Fi_j;
 end;



end.
