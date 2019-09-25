unit UnitUserDefinedFluidMaterial;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormUserDefinedFluidMaterial = class(TForm)
    PanelFLUID: TPanel;
    GBUserProperties: TGroupBox;
    BApply: TButton;
    GBRho: TGroupBox;
    CBRho: TComboBox;
    ERho: TEdit;
    LRho: TLabel;
    GBVolexpans: TGroupBox;
    EBeta_T: TEdit;
    LSIBeta_T: TLabel;
    GBdynamviscosity: TGroupBox;
    EMu: TEdit;
    LSIMu: TLabel;
    GBconductivity: TGroupBox;
    ELam: TEdit;
    LSILam: TLabel;
    GBheatcapacity: TGroupBox;
    ECp: TEdit;
    LSICp: TLabel;
    GBmatname: TGroupBox;
    EMatName: TEdit;
    GBtippattern: TGroupBox;
    CBTipPattern: TComboBox;
    CBMu: TComboBox;
    BEditMu: TButton;
    CBLam: TComboBox;
    ButtonLam: TButton;
    ButtonCp: TButton;
    CBCp: TComboBox;
    procedure BApplyClick(Sender: TObject);
    procedure CBTipPatternChange(Sender: TObject);
    procedure CBRhoChange(Sender: TObject);
    procedure CBMuChange(Sender: TObject);
    procedure BEditMuClick(Sender: TObject);
    procedure ButtonCpClick(Sender: TObject);
    procedure ButtonLamClick(Sender: TObject);
    procedure CBCpChange(Sender: TObject);
    procedure CBLamChange(Sender: TObject);
  private
    { Private declarations }
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  FormUserDefinedFluidMaterial: TFormUserDefinedFluidMaterial;

implementation
 uses
     VisualUnit, addBlockUnit, UnitnonNewtonianFluid, Unitusertempdepend;
{$R *.dfm}

 procedure TFormUserDefinedFluidMaterial.patchstring(var s : String);
var
   i : Integer;
begin
if (FormatSettings.DecimalSeparator='.') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]=',') then s[i]:='.';
   end;
end;

if (FormatSettings.DecimalSeparator=',') then
begin
   for i:=1 to length(s) do
   begin
      if (s[i]='.') then s[i]:=',';
   end;
end;
end;

// Задание FLUID материала пользователем
// при нажатии на кнопку Apply.
procedure TFormUserDefinedFluidMaterial.BApplyClick(Sender: TObject);
var
 s : String;
begin
   with Laplas.body[Laplas.itek] do
    begin
       Laplas.workmat[imatid].namemat:=EMatName.Text; // имя материала

       s:=Trim(ERho.Text);
       patchstring(s);
       ERho.Text:=s;
       s:=Trim(ECp.Text);
       patchstring(s);
       ECp.Text:=s;
       s:=Trim(ELam.Text);
       patchstring(s);
       ELam.Text:=s;

       Laplas.workmat[imatid].rho:=StrToFloat(ERho.Text); // плотность
       //Laplas.workmat[imatid].cp:=StrToFloat(ECp.Text); // теплоёмкость
       Laplas.workmat[imatid].n_cp:=1;
       SetLength(Laplas.workmat[imatid].temp_cp, Laplas.workmat[imatid].n_cp);
       SetLength(Laplas.workmat[imatid].arr_cp, Laplas.workmat[imatid].n_cp);
       Laplas.workmat[imatid].temp_cp[0]:=20.0;
       Laplas.workmat[imatid].arr_cp[0]:= StrToFloat(ECp.Text); // удельная теплоёмкость при постоянном давлении.
       //Laplas.workmat[imatid].lambda:=StrToFloat(ELam.Text); // теплопроводность
       Laplas.workmat[imatid].n_lam:=1;
       SetLength(Laplas.workmat[imatid].temp_lam, Laplas.workmat[imatid].n_lam);
       SetLength(Laplas.workmat[imatid].arr_lam, Laplas.workmat[imatid].n_lam);
       Laplas.workmat[imatid].temp_lam[0]:=20.0;
       Laplas.workmat[imatid].arr_lam[0]:= StrToFloat(ELam.Text); // теплопроводность

       if (CBMu.ItemIndex=0) then
       begin
          Laplas.workmat[imatid].ilawmu:=0;
          s:=Trim(EMu.Text);
          patchstring(s);
          EMu.Text:=s;

          Laplas.workmat[imatid].mu:=StrToFloat(EMu.Text); // динамическая вязкость
       end
       else
       begin
          // в качестве стартового значения вязкости берём среднее значение.
          Laplas.workmat[imatid].mu:=0.5*(Laplas.workmat[imatid].mumin+Laplas.workmat[imatid].mumax);
          if (Laplas.workmat[imatid].ilawmu=0) then Application.MessageBox('Change happens only when you click Edit','Attantion!');
          CBMu.ItemIndex:=0;
       end;

       s:=Trim(EBeta_T.Text);
       patchstring(s);
       EBeta_T.Text:=s;

       Laplas.workmat[imatid].beta_t:=StrToFloat(EBeta_T.Text);
       if (CBRho.ItemIndex=1) then Laplas.workmat[imatid].bBoussinesq:=1 else Laplas.workmat[imatid].bBoussinesq:=0;
       Laplas.workmat[imatid].blibmat:=0; // это материал определённый пользователем
       Laplas.workmat[imatid].ilibident:=0; // это материал определённый пользователем
    end;
end;

// Выбор шаблона материала.
procedure TFormUserDefinedFluidMaterial.CBTipPatternChange(
  Sender: TObject);
begin
    BEditMu.Visible:=false;
    EMu.Visible:=true;
    LSIMu.Visible:=true;
    CBMu.ItemIndex:=0; // const
    if (CBTipPattern.ItemIndex>0) then
    begin
      EMatName.Text:=Laplas.libfluid[CBTipPattern.ItemIndex-1].name;
      ERho.Text:=FloatToStr(Laplas.libfluid[CBTipPattern.ItemIndex-1].rho);
      ECp.Text:=FloatToStr(Laplas.libfluid[CBTipPattern.ItemIndex-1].cp);
      ELam.Text:=FloatToStr(Laplas.libfluid[CBTipPattern.ItemIndex-1].lam);
      EMu.Text:=FloatToStr(Laplas.libfluid[CBTipPattern.ItemIndex-1].mu);
      EBeta_T.Text:=FloatToStr(Laplas.libfluid[CBTipPattern.ItemIndex-1].beta_t);
   end
   else
   begin
       // выбрано no pattern
       EMatName.Text:='';
       ERho.Text:='';
       ECp.Text:='';
       ELam.Text:='';
       EMu.Text:='';
       EBeta_T.Text:='';
   end;
end;


// Включение или выключение модели Обербека-Буссинеска
procedure TFormUserDefinedFluidMaterial.CBRhoChange(Sender: TObject);
begin
   case CBRho.ItemIndex of
     0 : begin
            // const rho
            GBVolexpans.Visible:=false;
         end;
     1 : begin
            // Boussinesq Model
            GBVolexpans.Visible:=true;
            //ShowMessage('Please, define Vol. expansion coef. and later  Operating Temperature ');
            //Laplas.MainMemo.Lines.Add('Please, define Vol. expansion coef. and later Operating Temperature ');
         end;
   end;
end;

// Смена типа задания удельной теплоёмкости с константы на кусочно-линейную.
procedure TFormUserDefinedFluidMaterial.CBCpChange(Sender: TObject);
begin
   case CBCp.ItemIndex of
     0 : begin
            // Constant
            ButtonCp.Visible:=false;
            ECp.Visible:=true;
            LSICp.Visible:=true;
         end;
     1 : begin
            // Piecewise
            ButtonCp.Visible:=true;
            ECp.Visible:=false;
            LSICp.Visible:=false;
         end;
   end;
end;

// Смена типа задания теплопроводности с константы на кусочно-линейную.
procedure TFormUserDefinedFluidMaterial.CBLamChange(Sender: TObject);
begin
    case CBLam.ItemIndex of
       0 : begin
              // Constant
              ButtonLam.Visible:=false;
              ELam.Visible:=true;
              LSILam.Visible:=true;
           end;
       1 : begin
              // Piecewise
              ButtonLam.Visible:=true;
              ELam.Visible:=false;
              LSILam.Visible:=false;
           end;
    end;
end;

// Выбор модели для вязкости
procedure TFormUserDefinedFluidMaterial.CBMuChange(Sender: TObject);
begin
   case CBMu.ItemIndex of
   0 : begin
          // const
          BEditMu.Visible:=false;
          EMu.Visible:=true;
          LSIMu.Visible:=true;
       end;
   1 : begin
          // non-Newtonian
          BEditMu.Visible:=true;
          EMu.Visible:=false;
          LSIMu.Visible:=false;
       end;
   end;
end;

procedure TFormUserDefinedFluidMaterial.BEditMuClick(Sender: TObject);
begin
   // вызывает форму для задания закона неньютоновской жидкоссти

   FormnonNewtonFluid.ImageOstwald_de_Vel.Visible:=false;
   FormnonNewtonFluid.Imagecaisson.Visible:=false;
   FormnonNewtonFluid.ImagePrandtl.Visible:=false;
   FormnonNewtonFluid.ImageCarreau.Visible:=false;
   FormnonNewtonFluid.ImagePowell_Eyring.Visible:=false;
   FormnonNewtonFluid.ImageWilliamson.Visible:=false;

   if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilawmu=0) then
   begin
      // power-law fluid
      // Оствальд де Вель
      FormnonNewtonFluid.CBlaw.ItemIndex:=0;
      FormnonNewtonFluid.ImageOstwald_de_Vel.Visible:=true;
      FormnonNewtonFluid.LabelB.Visible:=false;
      FormnonNewtonFluid.EditB.Visible:=false;
      FormnonNewtonFluid.LabelC.Visible:=false;
      FormnonNewtonFluid.EditC.Visible:=false;
      FormnonNewtonFluid.Labeln.Visible:=true;
      FormnonNewtonFluid.Editn.Visible:=true;
   end
   else
   begin
      FormnonNewtonFluid.CBlaw.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilawmu-1;
      case FormnonNewtonFluid.CBlaw.ItemIndex of
         0 : begin
            // power-law fluid
            // Оствальд де Вель
            FormnonNewtonFluid.ImageOstwald_de_Vel.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=false;
            FormnonNewtonFluid.EditB.Visible:=false;
            FormnonNewtonFluid.LabelC.Visible:=false;
            FormnonNewtonFluid.EditC.Visible:=false;
            FormnonNewtonFluid.Labeln.Visible:=true;
            FormnonNewtonFluid.Editn.Visible:=true;
         end;
     1 : begin
            // Кессон
            FormnonNewtonFluid.Imagecaisson.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=true;
            FormnonNewtonFluid.EditB.Visible:=true;
            FormnonNewtonFluid.LabelC.Visible:=false;
            FormnonNewtonFluid.EditC.Visible:=false;
            FormnonNewtonFluid.Labeln.Visible:=false;
            FormnonNewtonFluid.Editn.Visible:=false;
         end;
     2 : begin
            // Прандтль
            FormnonNewtonFluid.ImagePrandtl.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=true;
            FormnonNewtonFluid.EditB.Visible:=true;
            FormnonNewtonFluid.LabelC.Visible:=false;
            FormnonNewtonFluid.EditC.Visible:=false;
            FormnonNewtonFluid.Labeln.Visible:=false;
            FormnonNewtonFluid.Editn.Visible:=false;
         end;
     3 : begin
            // Carreau
            FormnonNewtonFluid.ImageCarreau.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=true;
            FormnonNewtonFluid.EditB.Visible:=true;
            FormnonNewtonFluid.LabelC.Visible:=true;
            FormnonNewtonFluid.EditC.Visible:=true;
            FormnonNewtonFluid.Labeln.Visible:=true;
            FormnonNewtonFluid.Editn.Visible:=true;
         end;
     4 : begin
            // Пауэлл-Эйринг
            FormnonNewtonFluid.ImagePowell_Eyring.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=true;
            FormnonNewtonFluid.EditB.Visible:=true;
            FormnonNewtonFluid.LabelC.Visible:=true;
            FormnonNewtonFluid.EditC.Visible:=true;
            FormnonNewtonFluid.Labeln.Visible:=false;
            FormnonNewtonFluid.Editn.Visible:=false;
         end;
     5 : begin
            // Уильямсон
            FormnonNewtonFluid.ImageWilliamson.Visible:=true;
            FormnonNewtonFluid.LabelB.Visible:=true;
            FormnonNewtonFluid.EditB.Visible:=true;
            FormnonNewtonFluid.LabelC.Visible:=true;
            FormnonNewtonFluid.EditC.Visible:=true;
            FormnonNewtonFluid.Labeln.Visible:=false;
            FormnonNewtonFluid.Editn.Visible:=false;
         end;
      end;
   end;
   FormnonNewtonFluid.Editmin.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mumin);
   FormnonNewtonFluid.Editmax.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mumax);
   FormnonNewtonFluid.EditA.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Amu);
   FormnonNewtonFluid.EditB.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Bmu);
   FormnonNewtonFluid.EditC.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].Cmu);
   FormnonNewtonFluid.Editn.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].degreennmu);
   FormnonNewtonFluid.ShowModal; // Вызов формы для задания свойств
end;

// Пользовательская кусочно линейная температуро-зависимая теплоёмкость.
procedure TFormUserDefinedFluidMaterial.ButtonCpClick(Sender: TObject);
var
   i_4 : Integer;
begin
   Formusertempdepend.Caption:='Solid specific heat';
   Formusertempdepend.Label1.Caption:='The following curve specification consists of';
   Formusertempdepend.Label2.Caption:='a list of temperature/sp_heat pairs, which define a';
   Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
   Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
   Formusertempdepend.Label5.Caption:='solid specific heat units J/(kg*K).';
   Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0;  // Градусы Цельсия
   Formusertempdepend.Memopiecewiseproperties.Clear;
   for i_4 := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp-1 do
   begin
      Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[i_4])+' '+FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[i_4]));
   end;
   Formusertempdepend.identifier:=0;
   Formusertempdepend.ShowModal;
end;

procedure TFormUserDefinedFluidMaterial.ButtonLamClick(Sender: TObject);
var
 i_4 : Integer;
begin
   Formusertempdepend.Caption:='Solid conductivity';
   Formusertempdepend.Label1.Caption:='The following curve specification consists of a list';
   Formusertempdepend.Label2.Caption:='of temperature/conductivity pairs, which define a';
   Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
   Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
   Formusertempdepend.Label5.Caption:='solid conductivity units W/(m*K).';
   Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0; // Градусы Цельсия
   Formusertempdepend.Memopiecewiseproperties.Clear;
   for i_4 := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam-1 do
   begin
      Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[i_4])+' '+FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[i_4]));
   end;
   Formusertempdepend.identifier:=1;
   Formusertempdepend.ShowModal;
end;

end.
