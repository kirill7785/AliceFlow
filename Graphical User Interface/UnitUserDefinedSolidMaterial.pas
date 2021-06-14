unit UnitUserDefinedSolidMaterial;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormUserDefinedSolidMat = class(TForm)
    PanelSOLID: TPanel;
    GBuserproperties: TGroupBox;
    BApply: TButton;
    GroupBoxRho: TGroupBox;
    Erho: TEdit;
    LSIrho: TLabel;
    GroupBoxCp: TGroupBox;
    ECp: TEdit;
    LSICp: TLabel;
    GroupBoxLam: TGroupBox;
    ELam: TEdit;
    LSILam: TLabel;
    GroupBoxsmn: TGroupBox;
    EMatName: TEdit;
    GroupBoxtippattern: TGroupBox;
    CBsolidmat: TComboBox;
    GroupBoxOrthotropy: TGroupBox;
    Editmultx: TEdit;
    Editmulty: TEdit;
    Editmultz: TEdit;
    Labelx: TLabel;
    Labely: TLabel;
    Labelz: TLabel;
    ComboBoxconductivitytype: TComboBox;
    Buttonconductiviypiecewise: TButton;
    ComboBoxheatcapacitytype: TComboBox;
    Buttonheatcapacitypiecewise: TButton;
    GroupBoxThermalStress: TGroupBox;
    GroupBoxPoissonRatio: TGroupBox;
    GroupBoxYoungModule: TGroupBox;
    GroupBoxLinearExpansionCoefficient: TGroupBox;
    EditPoissonRatio: TEdit;
    EditYoungModule: TEdit;
    LabelYoungModule: TLabel;
    EditLinearExpansionKoefficient: TEdit;
    LabelLinearExpansionCoefficient: TLabel;
    ButtonCancel: TButton;
    ComboBoxlinearExpansion: TComboBox;
    ButtonEditlinearexpansioncoefficient: TButton;
    ButtonYoungModule: TButton;
    ComboBoxYoungModule: TComboBox;
    ComboBoxPoissonratio: TComboBox;
    ButtonPoissonratio: TButton;
    Label1: TLabel;
    EditbetaX: TEdit;
    Label2: TLabel;
    EditbetaY: TEdit;
    Label3: TLabel;
    EditbetaZ: TEdit;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Editnuxy: TEdit;
    Label7: TLabel;
    Editnuxz: TEdit;
    Label8: TLabel;
    Editnuyz: TEdit;
    GroupBoxShearModulus: TGroupBox;
    CheckBoxShearModulus: TCheckBox;
    PanelShearModulus: TPanel;
    Label9: TLabel;
    EditGxy: TEdit;
    Label10: TLabel;
    EditGyz: TEdit;
    Label11: TLabel;
    EditGxz: TEdit;
    Label12: TLabel;
    Label13: TLabel;
    EditEx: TEdit;
    Label14: TLabel;
    EditEy: TEdit;
    Label15: TLabel;
    EditEz: TEdit;
    procedure BApplyClick(Sender: TObject);
    procedure CBsolidmatChange(Sender: TObject);
    procedure ComboBoxheatcapacitytypeChange(Sender: TObject);
    procedure ComboBoxconductivitytypeChange(Sender: TObject);
    procedure ButtonheatcapacitypiecewiseClick(Sender: TObject);
    procedure ButtonconductiviypiecewiseClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ButtonCancelClick(Sender: TObject);
    procedure ButtonEditlinearexpansioncoefficientClick(Sender: TObject);
    procedure ComboBoxlinearExpansionChange(Sender: TObject);
    procedure ButtonYoungModuleClick(Sender: TObject);
    procedure ButtonPoissonratioClick(Sender: TObject);
    procedure ComboBoxYoungModuleClick(Sender: TObject);
    procedure ComboBoxPoissonratioClick(Sender: TObject);
    procedure CheckBoxShearModulusClick(Sender: TObject);

  private
    { Private declarations }

    procedure patchstring(var s : String);
  public
    { Public declarations }
    bCanselSelect : Boolean;

  end;

var
  FormUserDefinedSolidMat: TFormUserDefinedSolidMat;

implementation
uses
     VisualUnit, addBlockUnit, Unitusertempdepend;
{$R *.dfm}

procedure TFormUserDefinedSolidMat.patchstring(var s : String);
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


// Задаёт SOLID материал определённый пользователем
procedure TFormUserDefinedSolidMat.BApplyClick(Sender: TObject);
var
  s : String;
  bOk : Boolean;
begin

   bCanselSelect:=false;

    bOk:=true;
    // исправление десятичного сепаратора.
    s:=Trim(EditPoissonRatio.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       // Коэффициент Пуассона
       s:='0.154'; // AlSiC8
       patchstring(s);
    end;
    EditPoissonRatio.Text:=s;
    s:=Trim(EditYoungModule.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       // модуль Юнга.
       s:='217.5';// GPa AlSiC8
       patchstring(s);
    end;
    EditYoungModule.Text:=s;
    s:=Trim(EditLinearExpansionKoefficient.Text);
    if (length(Trim(s))=0) then
    begin
       // Коэффициент линейного теплового расширения.
       s:='6.5';// 1E-6 AlSiC8
       patchstring(s);
    end;
    patchstring(s);
    EditLinearExpansionKoefficient.Text:=s;

    s:=Trim(Erho.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='2800'; // density
       bOk:=false;
       patchstring(s);
    end;
    Erho.Text:=s;
    s:=Trim(ECp.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='921'; // heat capacity
       bOk:=false;
       patchstring(s);
    end;
    ECp.Text:=s;
    s:=Trim(ELam.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='164'; // thermal conductivity
       bOk:=false;
       patchstring(s);
    end;
    ELam.Text:=s;
    s:=Trim(Editmultx.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='1.0';
       bOk:=false;
       patchstring(s);
    end;
    Editmultx.Text:=s;
    s:=Trim(Editmulty.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='1.0';
       bOk:=false;
       patchstring(s);
    end;
    Editmulty.Text:=s;
    s:=Trim(Editmultz.Text);
    patchstring(s);
    if (length(Trim(s))=0) then
    begin
       s:='1.0';
       bOk:=false;
       patchstring(s);
    end;
    Editmultz.Text:=s;

    if (length(Trim(EMatName.Text))=0) then
    begin
        bOk:=false;
        EMatName.Text:='noname';
    end;

    if (bOk) then
    begin
       with Laplas.body[Laplas.itek] do
       begin
          Laplas.workmat[imatid].namemat:=EMatName.Text; // имя материала
          Laplas.workmat[imatid].rho:=StrToFloat(Erho.Text); // плотность
          //Laplas.workmat[imatid].cp:=StrToFloat(ECp.Text); // теплоёмкость
          if (ComboBoxHeatCapacitytype.ItemIndex=0) then
          begin
             // Constant properties
             Laplas.workmat[imatid].n_cp:=1;
             SetLength(Laplas.workmat[imatid].temp_cp, Laplas.workmat[imatid].n_cp);
             SetLength(Laplas.workmat[imatid].arr_cp, Laplas.workmat[imatid].n_cp);
             Laplas.workmat[imatid].temp_cp[0]:=20.0;
             Laplas.workmat[imatid].arr_cp[0]:= StrToFloat(ECp.Text); // удельная теплоёмкость при постоянном давлении
             // Если же свойства температурно зависимые при
             // ComboBoxHeatCapacitytype.ItemIndex=1
             // то уже все данные занесены и здесь ввод не нужен.
          end;
          //Laplas.workmat[imatid].lambda:=StrToFloat(ELam.Text); // теплопроводность
          if (ComboBoxconductivitytype.ItemIndex=0) then
          begin
             // Constant properties
             Laplas.workmat[imatid].n_lam:=1;
             SetLength(Laplas.workmat[imatid].temp_lam, Laplas.workmat[imatid].n_lam);
             SetLength(Laplas.workmat[imatid].arr_lam, Laplas.workmat[imatid].n_lam);
             Laplas.workmat[imatid].temp_lam[0]:=20.0;
             Laplas.workmat[imatid].arr_lam[0]:= StrToFloat(ELam.Text); // теплопроводность
             // Если же свойства температурно зависимые при
             // ComboBoxconductivitytype.ItemIndex=1
             // то уже все данные занесены и здесь ввод не нужен.
          end;
          // Ортотропность коэффициента теплопроводности.
          Laplas.workmat[imatid].mult_lam_x:=StrToFloat(Editmultx.Text);
          Laplas.workmat[imatid].mult_lam_y:=StrToFloat(Editmulty.Text);
          Laplas.workmat[imatid].mult_lam_z:=StrToFloat(Editmultz.Text);
          // Коэффициент линейного теплового расширения.
          Laplas.workmat[imatid].mult_Linear_expansion_coefficient_x:=StrToFloat(EditbetaX.Text);
          Laplas.workmat[imatid].mult_Linear_expansion_coefficient_y:=StrToFloat(EditbetaY.Text);
          Laplas.workmat[imatid].mult_Linear_expansion_coefficient_z:=StrToFloat(EditbetaZ.Text);
          // Модуль Юнга.
          Laplas.workmat[imatid].mult_Young_Module_x:=StrToFloat(EditEx.Text);
          Laplas.workmat[imatid].mult_Young_Module_y:=StrToFloat(EditEy.Text);
          Laplas.workmat[imatid].mult_Young_Module_z:=StrToFloat(EditEz.Text);
          // Ортотропность коэффициента Пуассона.
          Laplas.workmat[imatid].mult_Poisson_ratio_xy:=StrToFloat(Editnuxy.Text);
          Laplas.workmat[imatid].mult_Poisson_ratio_xz:=StrToFloat(Editnuxz.Text);
          Laplas.workmat[imatid].mult_Poisson_ratio_yz:=StrToFloat(Editnuyz.Text);
          // Для выполнения симметричности матрицы механических свойств относительно главной диагонали должны выполняться
          // следующие условия симметрии.
          Laplas.workmat[imatid].mult_Poisson_ratio_yx:=StrToFloat(Editnuxy.Text)*(Laplas.workmat[imatid].mult_Young_Module_y/Laplas.workmat[imatid].mult_Young_Module_x);
          Laplas.workmat[imatid].mult_Poisson_ratio_zx:=StrToFloat(Editnuxz.Text)*(Laplas.workmat[imatid].mult_Young_Module_z/Laplas.workmat[imatid].mult_Young_Module_x);
          Laplas.workmat[imatid].mult_Poisson_ratio_zy:=StrToFloat(Editnuyz.Text)*(Laplas.workmat[imatid].mult_Young_Module_z/Laplas.workmat[imatid].mult_Young_Module_y);
          // Модуль сдвига
          if (CheckBoxShearModulus.Checked) then
          begin
              Laplas.workmat[imatid].bShearModuleActive:=true;
          end
          else
          begin
             Laplas.workmat[imatid].bShearModuleActive:=false;
          end;
          Laplas.workmat[imatid].ShearModuleGxy:=StrToFloat(EditGxy.Text);
          Laplas.workmat[imatid].ShearModuleGyz:=StrToFloat(EditGyz.Text);
          Laplas.workmat[imatid].ShearModuleGxz:=StrToFloat(EditGxz.Text);


          // Thermal-Stress
          //Laplas.workmat[imatid].Poisson_ratio:= StrToFloat(EditPoissonRatio.Text);
           if (ComboBoxPoissonratio.ItemIndex=0) then
          begin
             // Constant properties
             Laplas.workmat[imatid].n_Poisson_ratio:=1;
             SetLength(Laplas.workmat[imatid].temp_Poisson_ratio, Laplas.workmat[imatid].n_Poisson_ratio);
             SetLength(Laplas.workmat[imatid].arr_Poisson_ratio, Laplas.workmat[imatid].n_Poisson_ratio);
             Laplas.workmat[imatid].temp_Poisson_ratio[0]:=20.0;
             Laplas.workmat[imatid].arr_Poisson_ratio[0]:= StrToFloat(EditPoissonRatio.Text); // Коэффициент Пуассона.
             // Если же свойства температурно зависимые при
             // ComboBoxPoissonratio.ItemIndex=1
             // то уже все данные занесены и здесь ввод не нужен.
          end;
          //Laplas.workmat[imatid].Young_Module:=StrToFloat(EditYoungModule.Text);
           if (ComboBoxYoungModule.ItemIndex=0) then
          begin
             // Constant properties
             Laplas.workmat[imatid].n_Young_Module:=1;
             SetLength(Laplas.workmat[imatid].temp_Young_Module, Laplas.workmat[imatid].n_Young_Module);
             SetLength(Laplas.workmat[imatid].arr_Young_Module, Laplas.workmat[imatid].n_Young_Module);
             Laplas.workmat[imatid].temp_Young_Module[0]:=20.0;
             Laplas.workmat[imatid].arr_Young_Module[0]:= StrToFloat(EditYoungModule.Text); // модуль Юнга.
             // Если же свойства температурно зависимые при
             // ComboBoxYoungModule.ItemIndex=1
             // то уже все данные занесены и здесь ввод не нужен.
          end;
          //Laplas.workmat[imatid].Linear_expansion_coefficient:=StrToFloat(EditLinearExpansionKoefficient.Text);
           if (ComboBoxlinearExpansion.ItemIndex=0) then
          begin
             // Constant properties
             Laplas.workmat[imatid].n_Linear_expansion_coefficient:=1;
             SetLength(Laplas.workmat[imatid].temp_Linear_expansion_coefficient, Laplas.workmat[imatid].n_Linear_expansion_coefficient);
             SetLength(Laplas.workmat[imatid].arr_Linear_expansion_coefficient, Laplas.workmat[imatid].n_Linear_expansion_coefficient);
             Laplas.workmat[imatid].temp_Linear_expansion_coefficient[0]:=20.0;
             Laplas.workmat[imatid].arr_Linear_expansion_coefficient[0]:= StrToFloat(EditLinearExpansionKoefficient.Text); // коэффициент линейного теплового расширения.
             // Если же свойства температурно зависимые при
             // ComboBoxlinearExpansion.ItemIndex=1
             // то уже все данные занесены и здесь ввод не нужен.
          end;
          Laplas.workmat[imatid].blibmat:=0; // это материал определённый пользователем
          Laplas.workmat[imatid].ilibident:=100; // это материал определённый пользователем
       end;
       Close;
    end;
end;

// Вызов piecewise thermal conductivity.
procedure TFormUserDefinedSolidMat.ButtonCancelClick(Sender: TObject);
begin
   bCanselSelect:=true;
   Close;// зактыть форму.
end;

procedure TFormUserDefinedSolidMat.ButtonconductiviypiecewiseClick(
  Sender: TObject);
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

procedure TFormUserDefinedSolidMat.ButtonEditlinearexpansioncoefficientClick(
  Sender: TObject);
  var
    i_4 : Integer;
begin
   Formusertempdepend.Caption:='Solid linear expansion coefficient';
   Formusertempdepend.Label1.Caption:='The following curve specification consists of';
   Formusertempdepend.Label2.Caption:='a list of temperature/linear expansion coefficient pairs, which define a';
   Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
   Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
   Formusertempdepend.Label5.Caption:='solid linear expansion coefficient units 1E-6/(K).';
   Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0;  // Градусы Цельсия
   Formusertempdepend.Memopiecewiseproperties.Clear;
   for i_4 := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Linear_expansion_coefficient-1 do
   begin
      Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_Linear_expansion_coefficient[i_4])+' '+FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Linear_expansion_coefficient[i_4]));
   end;
   Formusertempdepend.identifier:=3;
   Formusertempdepend.ShowModal;
end;

// Вызов piecewise heat capacity.
procedure TFormUserDefinedSolidMat.ButtonheatcapacitypiecewiseClick(
  Sender: TObject);
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

procedure TFormUserDefinedSolidMat.ButtonPoissonratioClick(Sender: TObject);
var
  i_4 : Integer;
begin
   Formusertempdepend.Caption:='Solid Poisson ratio ';
   Formusertempdepend.Label1.Caption:='The following curve specification consists of';
   Formusertempdepend.Label2.Caption:='a list of temperature/Poisson ratio pairs, which define a';
   Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
   Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
   Formusertempdepend.Label5.Caption:='solid Poisson ratio units 1.';
   Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0;  // Градусы Цельсия
   Formusertempdepend.Memopiecewiseproperties.Clear;
   for i_4 := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Poisson_ratio-1 do
   begin
      Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_Poisson_ratio[i_4])+' '+FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Poisson_ratio[i_4]));
   end;
   Formusertempdepend.identifier:=5;
   Formusertempdepend.ShowModal;
end;

procedure TFormUserDefinedSolidMat.ButtonYoungModuleClick(Sender: TObject);
var
   i_4 : Integer;
begin
   Formusertempdepend.Caption:='Solid Young Module GPa';
   Formusertempdepend.Label1.Caption:='The following curve specification consists of';
   Formusertempdepend.Label2.Caption:='a list of temperature/Young Module pairs, which define a';
   Formusertempdepend.Label3.Caption:='piecewise-linear curve. Spacing is not significant';
   Formusertempdepend.Label4.Caption:='as long as the numbers are given in pairs.';
   Formusertempdepend.Label5.Caption:='solid Young Module units GPa.';
   Formusertempdepend.ComboBoxtemperatureUnit.ItemIndex:=0;  // Градусы Цельсия
   Formusertempdepend.Memopiecewiseproperties.Clear;
   for i_4 := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_Young_Module-1 do
   begin
      Formusertempdepend.Memopiecewiseproperties.Lines.Add(FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_Young_Module[i_4])+' '+FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_Young_Module[i_4]));
   end;
   Formusertempdepend.identifier:=4;
   Formusertempdepend.ShowModal;
end;

// Выбор шаблона заполнения:
procedure TFormUserDefinedSolidMat.CBsolidmatChange(Sender: TObject);
begin
   if (CBsolidmat.ItemIndex>0) then
   begin
      EMatName.Text:=Laplas.libsolid[CBsolidmat.ItemIndex-1].name;
      Erho.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].rho);
      ECp.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].cp);
      ELam.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].lam);
      Editmultx.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].mult_lam_x);
      Editmulty.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].mult_lam_y);
      Editmultz.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].mult_lam_z);
      // Множитель коэффициента линейного теплового расширения
      EditbetaX.Text:='1.0';
      EditbetaY.Text:='1.0';
      EditbetaZ.Text:='1.0';
      // Множитель модуля Юнга.
      EditEx.Text:='1.0';
      EditEy.Text:='1.0';
      EditEz.Text:='1.0';
      // Множитель коэффициента Пуассона
      Editnuxy.Text:='1.0';
      Editnuxz.Text:='1.0';
      Editnuyz.Text:='1.0';

      CheckBoxShearModulus.Checked:=false;
      CheckBoxShearModulusClick(Sender);
      EditGxy.Text:='1.0';
      EditGyz.Text:='1.0';
      EditGxz.Text:='1.0';

      EditPoissonRatio.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].Poisson_Ratio);
      EditYoungModule.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].Young_Module);
      EditLinearExpansionKoefficient.Text:=FloatToStr(Laplas.libsolid[CBsolidmat.ItemIndex-1].Linear_expansion_coefficient);
      //ComboBoxconductivitytype.ItemIndex:=0;
      //ComboBoxheatcapacitytype.ItemIndex:=0;
      //ComboBoxheatcapacitytypeChange(Sender);
      //ComboBoxconductivitytypeChange(Sender);
   end
   else
   begin
       // выбрано no pattern
       EMatName.Text:='';
       Erho.Text:='';
       ECp.Text:='';
       ELam.Text:='';
       Editmultx.Text:='';
       Editmulty.Text:='';
       Editmultz.Text:='';
       EditPoissonRatio.Text:='';
       EditYoungModule.Text:='';
       EditLinearExpansionKoefficient.Text:='';
   end;
end;

procedure TFormUserDefinedSolidMat.CheckBoxShearModulusClick(Sender: TObject);
begin
   if (CheckBoxShearModulus.Checked) then
   begin
      PanelShearModulus.Visible:=true;
   end
   else
   begin
      PanelShearModulus.Visible:=false;
   end;
end;

// Смена типа задания теплопроводности с константы на кусочно-линейную.
procedure TFormUserDefinedSolidMat.ComboBoxconductivitytypeChange(
  Sender: TObject);
begin
    case ComboBoxconductivitytype.ItemIndex of
       0 : begin
              // Constant
              Buttonconductiviypiecewise.Visible:=false;
              ELam.Visible:=true;
              LSILam.Visible:=true;
           end;
       1 : begin
              // Piecewise
              Buttonconductiviypiecewise.Visible:=true;
              ELam.Visible:=false;
              LSILam.Visible:=false;
           end;
    end;
end;

// Смена типа задания удельной теплоёмкости с константы на кусочно-линейную.
procedure TFormUserDefinedSolidMat.ComboBoxheatcapacitytypeChange(
  Sender: TObject);
begin
   case ComboBoxheatcapacitytype.ItemIndex of
     0 : begin
            // Constant
            Buttonheatcapacitypiecewise.Visible:=false;
            ECp.Visible:=true;
            LSICp.Visible:=true;
         end;
     1 : begin
            // Piecewise
            Buttonheatcapacitypiecewise.Visible:=true;
            ECp.Visible:=false;
            LSICp.Visible:=false;
         end;
   end;
end;

procedure TFormUserDefinedSolidMat.ComboBoxlinearExpansionChange(
  Sender: TObject);
begin
   case ComboBoxlinearExpansion.ItemIndex of
       0 : begin
              // Constant
              ButtonEditlinearexpansioncoefficient.Visible:=false;
              EditLinearExpansionKoefficient.Visible:=true;
              LabelLinearExpansionCoefficient.Visible:=true;
           end;
       1 : begin
              // Piecewise
              ButtonEditlinearexpansioncoefficient.Visible:=true;
              EditLinearExpansionKoefficient.Visible:=false;
              LabelLinearExpansionCoefficient.Visible:=false;
           end;
    end;
end;

procedure TFormUserDefinedSolidMat.ComboBoxPoissonratioClick(Sender: TObject);
begin
   case ComboBoxPoissonratio.ItemIndex of
       0 : begin
              // Constant
              ButtonPoissonratio.Visible:=false;
              EditPoissonRatio.Visible:=true;
           end;
       1 : begin
              // Piecewise
              ButtonPoissonratio.Visible:=true;
              EditPoissonRatio.Visible:=false;
           end;
   end;
end;

procedure TFormUserDefinedSolidMat.ComboBoxYoungModuleClick(Sender: TObject);
begin
   case ComboBoxYoungModule.ItemIndex of
       0 : begin
              // Constant
              ButtonYoungModule.Visible:=false;
              EditYoungModule.Visible:=true;
              LabelYoungModule.Visible:=true;
           end;
       1 : begin
              // Piecewise
              ButtonYoungModule.Visible:=true;
              EditYoungModule.Visible:=false;
              LabelYoungModule.Visible:=false;
           end;
   end;
end;

procedure TFormUserDefinedSolidMat.FormCreate(Sender: TObject);
begin
    bCanselSelect:=false;
end;

end.
