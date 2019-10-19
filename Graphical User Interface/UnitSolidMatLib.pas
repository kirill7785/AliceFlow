unit UnitSolidMatLib;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, CheckLst;

type
  TFormSolidLibMat = class(TForm)
    GBSelectMat: TGroupBox;
    BApply: TButton;
    LMessage1: TLabel;
    cbbSolidMatLib: TComboBox;
    lblConduct: TLabel;
    lblcapacity: TLabel;
    lblDensity: TLabel;
    btnviewL: TButton;
    btncapacity: TButton;
    lblvalDensity: TLabel;
    lblDSI: TLabel;
    GroupBoxOrthotropy: TGroupBox;
    Editx: TEdit;
    Edity: TEdit;
    Editz: TEdit;
    Labelx: TLabel;
    Labely: TLabel;
    Labelz: TLabel;
    procedure BApplyClick(Sender: TObject);
    procedure btnviewLClick(Sender: TObject);
    procedure btncapacityClick(Sender: TObject);
    procedure cbbSolidMatLibChange(Sender: TObject);
  private
    { Private declarations }
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  FormSolidLibMat: TFormSolidLibMat;

implementation
uses
     VisualUnit, addBlockUnit, Unitrectangularplot;
{$R *.dfm}

procedure TFormSolidLibMat.patchstring(var s : String);
var
   i : Integer;
begin
if (FormatSettings.DecimalSeparator='.') then
begin
   (*
   for i:=1 to length(s) do
   begin
      if (s[i]=',') then s[i]:='.';
   end;
   *)
    s:=StringReplace(s,',','.',[rfReplaceAll]);
end;

if (FormatSettings.DecimalSeparator=',') then
begin
   (*
   for i:=1 to length(s) do
   begin
      if (s[i]='.') then s[i]:=',';
   end;
   *)
    s:=StringReplace(s,'.',',',[rfReplaceAll]);
end;
end;


// Нажатие на кнопку Apply-применить
// Выбор билиотечного материала.
procedure TFormSolidLibMat.BApplyClick(Sender: TObject);
var
   s : String;
   bOk : Boolean;
begin

  bOk:=true;

  // Выделен только один библиотечный материал
  with Laplas.body[Laplas.itek] do
  begin
     //
     Laplas.workmat[imatid].blibmat:=1; // означает что это библиотечный материал
     Laplas.workmat[imatid].ilibident:=101+cbbSolidMatLib.ItemIndex; // номер библиотечного материала в соответствии с заданным
     s:=Editx.Text;
     patchstring(s);
     if (length(Trim(s))=0) then
     begin
       s:='1.0';
       bOk:=false;
       patchstring(s);
     end;
     Editx.Text:=s;
     Laplas.workmat[imatid].mult_lam_x:=StrToFloat(Editx.Text);
     s:=Edity.Text;
     patchstring(s);
     if (length(Trim(s))=0) then
     begin
       s:='1.0';
       bOk:=false;
       patchstring(s);
     end;
     Edity.Text:=s;
     Laplas.workmat[imatid].mult_lam_y:=StrToFloat(Edity.Text);
     s:=Editz.Text;
     patchstring(s);
     if (length(Trim(s))=0) then
     begin
        s:='1.0';
        bOk:=false;
        patchstring(s);
     end;
     Editz.Text:=s;
     Laplas.workmat[imatid].mult_lam_z:=StrToFloat(Editz.Text);
     // внутри кода AliceFlowv0_06.
  end;

  // Закрытие формы после успешного выполнения операций.
  if (bOk) then
  begin
     Close;
  end;
end;

// показывает график теплопроводности.
procedure TFormSolidLibMat.btnviewLClick(Sender: TObject);
begin
   frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
   case cbbSolidMatLib.ItemIndex of
   0 : begin
          // Alumina (поликор).
          // Справка : В AuSn золота именно 80% и 20% олова.
          frmRectangularPlot.cht1.Title.Text[0]:='Alumina, Al2O3 polycrystalline aluminum oxide conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // По данным программы TDIM-Master.
          frmRectangularPlot.cht1.Series[0].AddXY(200,82,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(298,46,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,32.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,24.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,18.9,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,13,'',clRed);
          // Справка - теплопроводность золото-оловянного припоя.
          // Нельзя рассчитывать теплопроводность припоя на основе свойств двух его
          // известных компонент  и процентного состава. Т.к. получается соединение
          // с неупорядоченной структурой имеются сильно перемешанные и неупорядоченные
          // вставки олова в золото, нет дальнего порядка : всё это приводит к тому что
          // теплопроводность припоя хуже теплопроводности каждой из его компонент.
          // Теплопроводность в данном случае экспериментальная величина.
          frmRectangularPlot.ShowModal;
       end;
   1 : begin
          // Si
          frmRectangularPlot.cht1.Title.Text[0]:='Si conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,266,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,156,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,105,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,80,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,64,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,52,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,43,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   2 : begin
          // GaAs
          frmRectangularPlot.cht1.Title.Text[0]:='GaAs conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,64.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,41.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(325,38,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(350,35.1,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(375,32.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,30.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(450,26.7,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,23.8,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,19.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,16.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,14.3,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   3 : begin
          // GaN
          frmRectangularPlot.cht1.Title.Text[0]:='GaN conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,268,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,150,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(325,134,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(350,120,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(375,109,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,99,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(450,84,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,72,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,56,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,45,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,37,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
    4 : begin
          // SiC4H
          frmRectangularPlot.cht1.Title.Text[0]:='SiC4H conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(3,6.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(4,21.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(5,42.4,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(6,78.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(7,154.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(8,231.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(9,307.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(10,385.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(20,1619.1,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(30,3315.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(40,3893.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(50,3855,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(60,3315.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(70,2883.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(80,2528.88,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(90,2235.9,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(100,2004.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(120,1603.7,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(140,1295.28,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(160,1040.85,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(180,778.71,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(200,693.9,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(220,593.67,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(240,516.57,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(260,454.89,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,424.05,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(280,400.92,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,370.08,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(198.15,717.96,'',clRed);
         // frmRectangularPlot.cht1.Series[0].AddXY(298.15,390.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(398.15,253.83,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(498.15,181.78,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(598.15,138.41,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   5 : begin
          // Sapphire
          frmRectangularPlot.cht1.Title.Text[0]:='Sapphire conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(298,42,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,32,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,20,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,13,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   6 : begin
          // Diamond
          frmRectangularPlot.cht1.Title.Text[0]:='Diamond conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(300,1250,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(360,1100,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(420,1020,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(480,950,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   7 : begin
          // МД40   40% Cu, 60% Mo.
          frmRectangularPlot.cht1.Title.Text[0]:='MD40 conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(273,254.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,252.8,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,247.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,241.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,235.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,230.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,225.4,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(900,220.1,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1000,214.9,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1200,208.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1400,202.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1600,198.0,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   8 : begin
          // Au
          frmRectangularPlot.cht1.Title.Text[0]:='Au conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(2,200,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(3,270,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(4,340,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(5,430,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(6,500,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(7,570,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(8,640,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(9,720,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(10,790,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(15,795,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(20,800,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(60,439,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(70,413,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(80,420,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(90,380,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(100,359,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(110,343,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(120,343,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(130,342,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(140,342,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(150,341,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(160,340,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(170,339,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(180,338,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(190,337,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(200,336,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(210,335,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(220,333,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(230,330,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(240,327,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(250,324,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(260,322,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,319,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(280,318,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(293,315,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,314,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,311,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,304,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,298,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,284,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   9 : begin
          // SiO2
          // Программа SYMMIC даёт величину 1.27 Вт/(м-К) для
          // термически выращенного слоя SiO2 из газовой фазы.
          frmRectangularPlot.cht1.Title.Text[0]:='SiO2 thermally-grown conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(298,1.27,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(498,1.27,'',clRed);
          // Эти данные неизвестно откуда взяты :
          //frmRectangularPlot.cht1.Series[0].AddXY(298,1.6,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(373,1.7,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(773,2.1,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   10 : begin
          // Cu
          frmRectangularPlot.cht1.Title.Text[0]:='Cu conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,413,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,393,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,379,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   11 : begin
          // Kovar
          frmRectangularPlot.cht1.Title.Text[0]:='Kovar conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(273,14.1,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,14.7,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,15.6,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,17.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(973,19.3,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   12 : begin
          // Brass LS 59-1-L
          frmRectangularPlot.cht1.Title.Text[0]:='Brass LS 59-1-L conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          //frmRectangularPlot.cht1.Series[0].AddXY(273,108.784,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(300,108.784,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,105,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,116,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,128,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,141,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,156,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,170,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(900,183,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   13 : begin
          // Al-Duralumin Д16
          frmRectangularPlot.cht1.Title.Text[0]:='Al-Duralumin conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          //frmRectangularPlot.cht1.Series[0].AddXY(273,164,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(300,164,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(150,90,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(298,117.236,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,120,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,139.797,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(473,146.545,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,163.293,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   14 : begin
          // Aluminium Nitride
          frmRectangularPlot.cht1.Title.Text[0]:='Aluminium Nitride conductivity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='conductivity, W/(m*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,780,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,319,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,195,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,100,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1000,49,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   end;

end;

// показывает график теплоёмкости.
procedure TFormSolidLibMat.btncapacityClick(Sender: TObject);
begin
   frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;
   case cbbSolidMatLib.ItemIndex of
   0 : begin
          // Справка : AuSn   80% золота и 20% олова.
          frmRectangularPlot.cht1.Title.Text[0]:='Alumina Al2O3 polycrystalline aluminum oxide capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          // Программа TDIM - Master даёт значение теплоёмкости 143.
          // Таковы данные программы SYMMIC.
          frmRectangularPlot.cht1.Series[0].AddXY(200,502,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(298,753,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,920,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,1046,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,1088,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,1172,'',clRed);
           // Справка : Нельзя рассчитывать теплоёмкость припоя на основе свойств двух его
           // известных компонент  и процентного состава. Т.к. получается соединение
           // с неупорядоченной структурой имеются сильно перемешанные и неупорядоченные
           // вставки олова в золото, нет дальнего порядка : всё это приводит к тому что
           // теплоёмкость припоя хуже теплоёмкости каждой из его компонент.
           // Теплоёмкость в данном случае экспериментальная величина.
          frmRectangularPlot.ShowModal;
       end;
   1 : begin
          // Si
          frmRectangularPlot.cht1.Title.Text[0]:='Si capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(77,180,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(173,490,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,680,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,770,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,850,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,880,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   2 : begin
          // GaAs
          frmRectangularPlot.cht1.Title.Text[0]:='GaAs capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          //frmRectangularPlot.cht1.Series[0].AddXY(273,325,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(300,325,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(173,307.4,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(298.15,322.92,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(400,335.44,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(500,344.01,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(600,351.13,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(700,357.5,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(800,363.5,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(900,369.24,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1000,374.9,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1100,380.44,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1200,385.9,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1300,391.37,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1400,396.76,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1500,402.15,'',clRed);
           frmRectangularPlot.cht1.Series[0].AddXY(1514,402.84,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   3 : begin
          // GaN
          frmRectangularPlot.cht1.Title.Text[0]:='GaN capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          //frmRectangularPlot.cht1.Series[0].AddXY(273,490,'',clRed);
          //frmRectangularPlot.cht1.Series[0].AddXY(300,490,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(200,322.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,431.3,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,501.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,543.8,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,572.17,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(700,592.8,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,608.9,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(900,622.2,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1000,633.7,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1100,643.99,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1200,653.4,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1300,662.22,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
    4 : begin
          // SiC4H
          frmRectangularPlot.cht1.Title.Text[0]:='SiC4H capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(77,50,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,670,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,820,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,1010,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,1120,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   5 : begin
          // Sapphire
          frmRectangularPlot.cht1.Title.Text[0]:='Sapphire capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(77,60,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(173,403,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,718,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,907,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,1089,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,1168,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   6 : begin
          // Diamond
          frmRectangularPlot.cht1.Title.Text[0]:='Diamond capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(77,8,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(173,140,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,420,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,770,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,1300,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,1590,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   7 : begin
          // МД40     40% Cu 60% Mo
          frmRectangularPlot.cht1.Title.Text[0]:='MD40 capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(173,268.77,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,302.5,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,318.17,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,335.11,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   8 : begin
          // Au
          frmRectangularPlot.cht1.Title.Text[0]:='Au capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,124,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,129,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,131,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,133,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,135,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,140,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   9 : begin
          // SiO2
          // По данным программы SYMMIC данная зависимость теплоёмкости
          // от температуры справедлива для термически выращенного из
          // газовой фазы слоя SiO2.
          frmRectangularPlot.cht1.Title.Text[0]:='SiO2 capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(273,700,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,830,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,1020,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(773,1110,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   10 : begin
          // Cu
          frmRectangularPlot.cht1.Title.Text[0]:='Cu capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,356,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(300,385,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(400,397,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(500,412,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(600,417,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(800,433,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   11 : begin
          // Kovar
          frmRectangularPlot.cht1.Title.Text[0]:='Kovar capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(273.15,439.614,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(703.15,648.954,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   12 : begin
          // Brass LS 59-1-L
          frmRectangularPlot.cht1.Title.Text[0]:='Brass LS 59-1-L capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(77,200,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(173,340,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(273,387,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(373,390,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,448,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   13 : begin
          // Al-Duralumin Д16.
          frmRectangularPlot.cht1.Title.Text[0]:='Al-Duralumin capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(373,921,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(473,1047,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(573,1130,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(673,1172,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   14 : begin
          // Нитрид алюминия
          frmRectangularPlot.cht1.Title.Text[0]:='Aluminium Nitride capacity';
          frmRectangularPlot.cht1.LeftAxis.Title.Caption:='capacity, J/(kg*K)';
          frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, K';
          frmRectangularPlot.cht1.Series[0].Clear;
          frmRectangularPlot.cht1.Series[0].AddXY(200,600,'',clRed);
          frmRectangularPlot.cht1.Series[0].AddXY(1000,600,'',clRed);
          frmRectangularPlot.ShowModal;
       end;
   end;

end;

// На метке показывает плотность материала.
procedure TFormSolidLibMat.cbbSolidMatLibChange(Sender: TObject);
begin
   lblvalDensity.Caption:='';
   case cbbSolidMatLib.ItemIndex of
   0 : begin
          // Alumina Al2O3 polycrystalline aluminum oxide.
          lblvalDensity.Caption:='3960';
       end;
   1 : begin
          // Si
          lblvalDensity.Caption:='2330';
       end;
   2 : begin
          // GaAs
          lblvalDensity.Caption:='5316';
       end;
   3 : begin
          // GaN
          lblvalDensity.Caption:='6150';
       end;
   4 : begin
          // SiC4H
          lblvalDensity.Caption:='3200';
       end;
   5 : begin
          // Sapphire
          lblvalDensity.Caption:='3970';
       end;
   6 : begin
          // Diamond
          lblvalDensity.Caption:='3515';
       end;
   7 : begin
          // МД40
          lblvalDensity.Caption:='9600';
       end;
   8 : begin
          // Au
          lblvalDensity.Caption:='19300';
       end;
   9 : begin
          // SiO2
          lblvalDensity.Caption:='2100';
       end;
   10 : begin
          // Cu
          lblvalDensity.Caption:='8933';
       end;
   11 : begin
          // Kovar
          lblvalDensity.Caption:='8360';
       end;
   12 : begin
          // Латунь 59-1-Л
          lblvalDensity.Caption:='8500';
       end;
   13 : begin
          // Al-Duralumin
          lblvalDensity.Caption:='2800';
       end;
   14 : begin
          // Нитрид алюминия
          lblvalDensity.Caption:='3255';
       end;
   end;
   // Новые настройки всегда изотропны.
    FormSolidLibMat.Editx.Text:='1';
    FormSolidLibMat.Edity.Text:='1';
    FormSolidLibMat.Editz.Text:='1';

end;


end.
