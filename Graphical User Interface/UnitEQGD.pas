unit UnitEQGD;
// Аналог формы Basic Parameters в ANSYS Icepak 12.0.
// Логика настроек сильно упрощена
// 6 августа 2016 года.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, ExtCtrls, StdCtrls, Vcl.AppEvnts;

type
  TEGDForm = class(TForm)
    PanelGlobal: TPanel;
    GBCFlZ: TGroupBox;
    CBFlow: TCheckBox;
    RadioGroupFlowRegime: TRadioGroup;
    GroupBoxTurbulentModel: TGroupBox;
    ComboBoxturbulentmodel: TComboBox;
    ButtonidFlow: TButton;
    BEditTurb: TButton;
    GBGravity: TGroupBox;
    Lgx: TLabel;
    Lgy: TLabel;
    Lgz: TLabel;
    lblgx: TLabel;
    lblgy: TLabel;
    lbldz: TLabel;
    Egx: TEdit;
    Egy: TEdit;
    Egz: TEdit;
    ApplicationEvents1: TApplicationEvents;
    CheckBoxStaticStructural: TCheckBox;
    GroupBoxTemperature: TGroupBox;
    ComboBoxTemperature: TComboBox;
    //procedure ButtonTempClick(Sender: TObject);
    procedure RadioGroupFlowRegimeClick(Sender: TObject);
    procedure CBFlowClick(Sender: TObject);
    procedure ButtonidFlowClick(Sender: TObject);
    procedure CBIdCurFLzoneClick(Sender: TObject);
    procedure ComboBoxturbulentmodelChange(Sender: TObject);
    procedure BEditTurbClick(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure ComboBoxTemperatureChange(Sender: TObject);



  private
    { Private declarations }
    procedure patchstring(var s : String);
  public
    { Public declarations }
  end;

var
  EGDForm: TEGDForm;

implementation
uses
     VisualUnit, UnitSmagorinsky;

{$R *.dfm}

// Считывание информации о том нужно ли решать уравнение теплопроводности.
{*
procedure TEGDForm.ButtonTempClick(Sender: TObject);
begin
   Laplas.egddata.itemper:=ComboBoxTemperature.ItemIndex;
end;
*}

{*
// Реакция на смену глобального числа FLUID зон.
// Данный код устарел т.к. 5 августа 2016 все гидродинамические подобласти
// слиты в одну и решатель при этом считал.
procedure TEGDForm.BmaxFluidDomainClick(Sender: TObject);
var
   i : Integer;
begin
   // считываем количество максимальное FLUID зон
   Laplas.egddata.imaxflD:=CBMaxFluidDomain.ItemIndex;
   // Выдиление или уничтожение оперативной памяти из-за изменения размеров
   // динамического массива.
   SetLength(Laplas.egddata.myflmod,CBMaxFluidDomain.ItemIndex);

   // инициализация нулём
   for i:=0 to Laplas.egddata.imaxflD-1 do
   begin
      Laplas.egddata.myflmod[i].xc:=0.0;
      Laplas.egddata.myflmod[i].yc:=0.0;
      Laplas.egddata.myflmod[i].zc:=0.0;
      Laplas.egddata.myflmod[i].iflow:=0; // не считаем течение
      Laplas.egddata.myflmod[i].iflowregime:=0; // laminar
      Laplas.egddata.myflmod[i].iturbmodel:=0; // 0-ZEM, 1-Smagorinsky, 2 - RNG LES.
      // модель Смагоринского
      BEditTurb.Visible:=false; // делаем кнопку редатирования свойств модели Смагоринского невидимой, т.к. активна ZEM
      Laplas.egddata.myflmod[i].SmagConst:=0.151; // при Ck==1.8   (Ck соответствующая константа Колмогорова).
      Laplas.egddata.myflmod[i].bSmagorinsky_Lilly:=true; // модель Смагоринского-Лиллу.
      Laplas.egddata.myflmod[i].bfdelta:=true; // учёт неравномерности сетки
      Laplas.egddata.myflmod[i].bsurface_roughness:=false; // не учитывать шероховатость стенки.
      Laplas.egddata.myflmod[i].ipowerroughness:=2; // показатель степени в модели.
      Laplas.egddata.myflmod[i].roughness:=10.0; // micron (шероховатость стенки в мкм).
      Laplas.egddata.myflmod[i].bSwirlamendment:=false;
      Laplas.egddata.myflmod[i].bSelectiveSmagorinsky:=False; // избирательная модель Смагоринского.
      //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr:=2; // фильтр Симпсона.
      //Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle:=15.0;
       Laplas.egddata.myflmod[0].itypefiltr:=2; // фильтр Симпсона.
       Laplas.egddata.myflmod[0].rSelectiveAngle:=15.0;
   end;
   Application.MessageBox('information on the hydrodynamic DOMAIN cleared','Attantion!');

   if (CBMaxFluidDomain.ItemIndex=0) then
   begin
      GBCFlZ.Visible:=false;
   end
    else
   begin
      GBCFlZ.Visible:=true;

     // CBIdCurFLzone.Clear;
      //for i:=0 to Laplas.egddata.imaxflD-1 do
      //begin
        // CBIdCurFLzone.Items.Append(IntToStr(i+1));
      //end;
      //CBIdCurFLzone.ItemIndex:=0; // устанавливаем первую FLUID зону

      // заполнение координат опорной точки
      //EditXC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc);
      //EditYC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc);
      //EditZC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc);

      // расчитывать или не расчитывать уравнения гидродинамики
     // if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow=1) then
      if (Laplas.egddata.myflmod[0].iflow=1) then
      begin
         CBFlow.Checked:=true;
      end
       else
      begin
         CBFlow.Checked:=false;
      end;

      // режим течения
      //RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime;
      // модель турбулентности
      //ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel;
      // режим течения
      RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[0].iflowregime;
      // модель турбулентности
      ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[0].iturbmodel;

   end;
end;
*}


procedure TEGDForm.RadioGroupFlowRegimeClick(Sender: TObject);
begin
   // смена режима течения:
   // ламинарный или турбулентный.
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      // LAMINAR
      GroupBoxTurbulentModel.Visible:=false;
   end
    else
   begin
      // Turbulent
      GroupBoxTurbulentModel.Visible:=true;
   end;
end;

procedure TEGDForm.CBFlowClick(Sender: TObject);
begin
   // решать или не решать уравнения гидродинамики
   if (CBFlow.Checked) then
   begin
      // уравнения гидродинамики решаются.
      if ((ComboBoxTemperature.ItemIndex=2) or
       (ComboBoxTemperature.ItemIndex=3)) then
      begin
         ComboBoxTemperature.ItemIndex:=1;
      end;

      RadioGroupFlowRegime.Visible:=true;
      GBGravity.Visible:=true;
      // ламинарный или турбулентный.
      if (RadioGroupFlowRegime.ItemIndex=0) then
      begin
         // LAMINAR
         GroupBoxTurbulentModel.Visible:=false;
      end
       else
      begin
         // Turbulent
         GroupBoxTurbulentModel.Visible:=true;
      end;
   end
    else
   begin
      // уравнения гидродинамики не решаются.
      RadioGroupFlowRegime.Visible:=false;
      GroupBoxTurbulentModel.Visible:=false;
      GBGravity.Visible:=false;
      if ((ComboBoxTemperature.ItemIndex=0) and
       (CheckBoxStaticStructural.Checked=false)) then
      begin
         // Не рассчитываем гидродинамику
         // Не рассчитываем механику.
         // Рассчитываем температуру методом
         // контрольного объёма.
         ComboBoxTemperature.ItemIndex:=1;
      end;

   end;
   // информация о том рассчитывать уравнения гидродинамики или нет.
   {*
   if (CBFlow.Checked) then
   begin
      // считаются уравнения гидродинамики.
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=1;
   end
    else
   begin
      // не считаются уравнения гидродинамики.
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=0;
   end;
   *}
    if (CBFlow.Checked) then
   begin
      // считаются уравнения гидродинамики.
      Laplas.egddata.myflmod[0].iflow:=1;
      // Расширяем форму чтобы показать гидродинамические
      // настройки решателя.
      Height:=488;
      ButtonidFlow.Top:=421;
      GBCFlZ.Height:=449;
      PanelGlobal.Height:=449;
   end
    else
   begin
      // не считаются уравнения гидродинамики.
      Laplas.egddata.myflmod[0].iflow:=0;
      // Скрываем ненужную в данный момент информацию
      // о гидродинамических характеристиках и
      // ужимаем форму для показа.
      GBCFlZ.Height:=135;
      Height:=136;
      PanelGlobal.Height:=136;
      ButtonidFlow.Top:=105;
   end;
end;

procedure TEGDForm.patchstring(var s : String);
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

procedure TEGDForm.ButtonidFlowClick(Sender: TObject);
var
  s1 : String;
  k : Integer;
   s : String;
  bOk : Boolean;
  c : Real;
  code : Integer;
  sforval : String;
begin
   // Начало считывания вектора силы тяжести.

   // исправление десятичного разделителя.
   s:=Trim(Egx.Text);
   patchstring(s);
   Egx.Text:=s;
   s:=Trim(Egy.Text);
   patchstring(s);
   Egy.Text:=s;
   s:=Trim(Egz.Text);
   patchstring(s);
   Egz.Text:=s;


   bOk:=False;
   //val(Trim(Egx.Text),c,code);
   sforval:='';
  sforval:=StringReplace(Egx.Text,',','.',[rfReplaceAll]);
   val(sforval,c,code);
   if (code=0) then
   begin
       //val(Trim(Egy.Text),c,code);
       sforval:='';
       sforval:=StringReplace(Egy.Text,',','.',[rfReplaceAll]);
       val(sforval,c,code);
      if (code=0) then
      begin
         //val(Trim(Egz.Text),c,code);
         sforval:='';
         sforval:=StringReplace(Egz.Text,',','.',[rfReplaceAll]);
         val(sforval,c,code);
         if (code=0) then
         begin
            bOk:=True;
         end;
      end;
   end;
   if (bOk) then
   begin
      Laplas.gx:=StrToFloat(Trim(Egx.Text));
      Laplas.gy:=StrToFloat(Trim(Egy.Text));
      Laplas.gz:=StrToFloat(Trim(Egz.Text));
   end
    else
   begin
      if (FormatSettings.DecimalSeparator='.') then
      begin
         Egx.Text:='0.0';
         Egy.Text:='0.0';
         Egz.Text:='0.0';
      end;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         Egx.Text:='0,0';
         Egy.Text:='0,0';
         Egz.Text:='0,0';
      end;
      ShowMessage('Input Error, please repeat...');
   end;

    // Конец считывания вектора силы тяжести.


   // считывание информации о текущей FLUID зоне
   {*
   s1:=Trim(EditXC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditXC.Text:=s1;

   s1:=Trim(EditYC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditYC.Text:=s1;

   s1:=Trim(EditZC.Text);
   for k:=1 to length(s1) do
   begin
      if (FormatSettings.DecimalSeparator=',') then
      begin
         if (s1[k]='.') then s1[k]:=',';
      end;
       if (FormatSettings.DecimalSeparator='.') then
      begin
         if (s1[k]=',') then s1[k]:='.';
      end;
   end;
   EditZC.Text:=s1;


   // Координаты опорной точки
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc:=StrToFloat(EditXC.Text);
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc:=StrToFloat(EditYC.Text);
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc:=StrToFloat(EditZC.Text);
    *}

    // Настройки для решения уравнения теплопередачи.
    // Считывание информации о том нужно ли решать уравнение теплопроводности.
    Laplas.egddata.itemper:=ComboBoxTemperature.ItemIndex;

   // Настройки для решения уравнений теории упругости.
   // Считывание информации о том нужно ли решать сиситему уравнений упругости.
   if (CheckBoxStaticStructural.Checked) then
   begin
      Laplas.egddata.iStaticStructural:=1;
   end
    else
   begin
      Laplas.egddata.iStaticStructural:=0;
   end;


    // С 5 августа 2016 года у нас тольк одна глобальная
    // FLUID зона, что позволяет считать многосвязные
    // гидродинамические подобласти.
    // Координаты опорной точки
   Laplas.egddata.myflmod[0].xc:=0.0;
   Laplas.egddata.myflmod[0].yc:=0.0;
   Laplas.egddata.myflmod[0].zc:=0.0;

   {*
   // информация о том рассчитывать уравнения гидродинамики или нет.
   if (CBFlow.Checked) then
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=1;
   end
    else
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow:=0;
   end;

   // режим течения 0-Ламинарный или 1-турбулентный
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime:=0;
   end
    else
   begin
      Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime:=1;
   end;

   // считывание информации о модели турбулентности.
   Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel:=ComboBoxturbulentmodel.ItemIndex;
   *}
   // информация о том рассчитывать уравнения гидродинамики или нет.
   if (CBFlow.Checked) then
   begin
      Laplas.egddata.myflmod[0].iflow:=1;
   end
    else
   begin
      Laplas.egddata.myflmod[0].iflow:=0;
   end;

   // режим течения 0-Ламинарный или 1-турбулентный
   if (RadioGroupFlowRegime.ItemIndex=0) then
   begin
      Laplas.egddata.myflmod[0].iflowregime:=0;
   end
    else
   begin
      Laplas.egddata.myflmod[0].iflowregime:=1;
   end;

   // считывание информации о модели турбулентности.
   Laplas.egddata.myflmod[0].iturbmodel:=ComboBoxturbulentmodel.ItemIndex;

   Close();

end;

procedure TEGDForm.CBIdCurFLzoneClick(Sender: TObject);
begin
   // смена текущей FLUID зоны

   // заполнение координат опорной точки
  // EditXC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].xc);
   //EditYC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].yc);
   //EditZC.Text:= FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].zc);

   // расчитывать или не расчитывать уравнения гидродинамики
   //if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflow=1) then
   if (Laplas.egddata.myflmod[0].iflow=1) then
   begin
      CBFlow.Checked:=true;
   end
    else
   begin
      CBFlow.Checked:=false;
   end;

   // режим течения
   //RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iflowregime;
   // модель турбулентности
   //ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iturbmodel;
    // режим течения
   RadioGroupFlowRegime.ItemIndex:=Laplas.egddata.myflmod[0].iflowregime;
   // модель турбулентности
   ComboBoxturbulentmodel.ItemIndex:=Laplas.egddata.myflmod[0].iturbmodel;
end;

// Делает видимой или невидимой кнопку редактирования
// Свойств модели Смагоринского.
procedure TEGDForm.ComboBoxTemperatureChange(Sender: TObject);
begin
   if ((ComboBoxTemperature.ItemIndex=0)
   and(CBFlow.Checked=false)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
       // Ставим по умолчанию солвер методом контрольного объёма.
       ComboBoxTemperature.ItemIndex:=1; // Не даём убирать температуру,
       // если выключена гидродинамика и механика.
   end;
   if ((ComboBoxTemperature.ItemIndex=2)and
   (CBFlow.Checked=true)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
      // Finite Element Method -> Finite Volume Method.
      ComboBoxTemperature.ItemIndex:=1;
   end;
    if ((ComboBoxTemperature.ItemIndex=3)and
   (CBFlow.Checked=true)and
   (CheckBoxStaticStructural.Checked=false)) then
   begin
      // Network T Solver -> Finite Volume Method.
      ComboBoxTemperature.ItemIndex:=1;
   end;
end;

procedure TEGDForm.ComboBoxturbulentmodelChange(Sender: TObject);
begin
   // Смена модели турбулентности
   case ComboBoxturbulentmodel.ItemIndex of
     0 : begin
            // ZEM (RANS)
            BEditTurb.Visible:=false;
         end;
     1 : begin
            // Smagorinsky
            BEditTurb.Visible:=true;
         end;
     2 : begin
            // RNG LES
            BEditTurb.Visible:=false;
         end;
     3 : begin
            // Spalart Allmares (RANS) [1992]
            BEditTurb.Visible:=false;
         end;
     4 : begin
            // SST Ментера (RANS) [1993]
            BEditTurb.Visible:=false;
         end;
     5 : begin
            // Двухслойная модель на основе
            // Standart K-Epsilon (RANS) [2001]
            BEditTurb.Visible:=false;
         end;
   end;
end;




// Инициализация формы GravityForm при создании формы.
procedure TEGDForm.FormCreate(Sender: TObject);
begin
    Egx.Text:=FloatToStr(Laplas.gx);
    Egy.Text:=FloatToStr(Laplas.gy);
    Egz.Text:=FloatToStr(Laplas.gz);
end;

// Запрет форме сворачиваться.
procedure TEGDForm.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

procedure TEGDForm.BEditTurbClick(Sender: TObject);
begin
{*
   // Изменение вызвано тем что 5 августа 2016 мы слили все гидродинамические подобласти
   // в одну единственную  и программа при этом работала.
   // Нажатие по кнопке редактирования модели Смагоринского.
   if (ComboBoxturbulentmodel.ItemIndex=1) then
   begin
      if (Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].iDynamicStressGermano=1) then
      begin
         // модель Германо включена.
         FormSmagorinsky.PanelCs.Visible:=false;
         FormSmagorinsky.PanelLimitersCs.Visible:=true;
         FormSmagorinsky.CBDynamicStress.Checked:=true;
         if (Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].iLimitersCs=0) then
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=false;
            FormSmagorinsky.Panel_user_limiters.Visible:=false;
         end
          else
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=true;
            // пределы Cs включены
            FormSmagorinsky.Panel_user_limiters.Visible:=true;
            FormSmagorinsky.EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].minCs);
            FormSmagorinsky.EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].maxCs);
         end;
      end
       else
      begin
         // модель Германо выключена.
         FormSmagorinsky.PanelCs.Visible:=true;
         FormSmagorinsky.Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].SmagConst);
         FormSmagorinsky.PanelLimitersCs.Visible:=false;
         FormSmagorinsky.CBDynamicStress.Checked:=false;
      end;

      FormSmagorinsky.CBcorrectnonuniformgrid.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bfdelta;
      FormSmagorinsky.CBSmagLilly.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bSmagorinsky_Lilly;
      FormSmagorinsky.CBsurfaceroughness.Checked:=Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].bsurface_roughness;
      if (FormSmagorinsky.CBsurfaceroughness.Checked) then
      begin
         // нужно активировать форму для задания шероховатости стенки.
         FormSmagorinsky.Panrougness.Visible:=true;
         FormSmagorinsky.Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].roughness);
         case Laplas.egddata.myflmod[CBIdCurFLzone.ItemIndex].ipowerroughness of
           1 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=0;
               end;
           2 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=1;
               end;
         end;
      end
      else
      begin
         FormSmagorinsky.Panrougness.Visible:=false;
      end;
      // поправка для течений с кривизной линий тока.
      FormSmagorinsky.CBSwirlAmendment.Checked:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSwirlamendment;
      if (FormSmagorinsky.CBSwirlAmendment.Checked) then
      begin
         FormSmagorinsky.PanRichardson.Visible:=true;
         FormSmagorinsky.ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rRimult);
      end
        else
      begin
         FormSmagorinsky.PanRichardson.Visible:=false;
      end;
      // Установка модели Selective Smagorinsky
      FormSmagorinsky.CBSelectiveSmagorinsky.Checked:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].bSelectiveSmagorinsky;
      if (FormSmagorinsky.CBSelectiveSmagorinsky.Checked) then
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=true;
         FormSmagorinsky.RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].itypefiltr;
         FormSmagorinsky.Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[EGDForm.CBIdCurFLzone.ItemIndex].rSelectiveAngle);
      end
       else
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=false;
      end;
      FormSmagorinsky.ShowModal;
   end;
   *}
   // Нажатие по кнопке редактирования модели Смагоринского.
   if (ComboBoxturbulentmodel.ItemIndex=1) then
   begin
      if (Laplas.egddata.myflmod[0].iDynamicStressGermano=1) then
      begin
         // модель Германо включена.
         FormSmagorinsky.PanelCs.Visible:=false;
         FormSmagorinsky.PanelLimitersCs.Visible:=true;
         FormSmagorinsky.CBDynamicStress.Checked:=true;
         if (Laplas.egddata.myflmod[0].iLimitersCs=0) then
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=false;
            FormSmagorinsky.Panel_user_limiters.Visible:=false;
         end
          else
         begin
            FormSmagorinsky.CBLimitersCS.Checked:=true;
            // пределы Cs включены
            FormSmagorinsky.Panel_user_limiters.Visible:=true;
            FormSmagorinsky.EminCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].minCs);
            FormSmagorinsky.EmaxCs.Text:=FloatToStr(Laplas.egddata.myflmod[0].maxCs);
         end;
      end
       else
      begin
         // модель Германо выключена.
         FormSmagorinsky.PanelCs.Visible:=true;
         FormSmagorinsky.Esmagconst.Text:=FloatToStr(Laplas.egddata.myflmod[0].SmagConst);
         FormSmagorinsky.PanelLimitersCs.Visible:=false;
         FormSmagorinsky.CBDynamicStress.Checked:=false;
      end;

      FormSmagorinsky.CBcorrectnonuniformgrid.Checked:=Laplas.egddata.myflmod[0].bfdelta;
      FormSmagorinsky.CBSmagLilly.Checked:=Laplas.egddata.myflmod[0].bSmagorinsky_Lilly;
      FormSmagorinsky.CBsurfaceroughness.Checked:=Laplas.egddata.myflmod[0].bsurface_roughness;
      if (FormSmagorinsky.CBsurfaceroughness.Checked) then
      begin
         // нужно активировать форму для задания шероховатости стенки.
         FormSmagorinsky.Panrougness.Visible:=true;
         FormSmagorinsky.Eroughnessmicron.Text:=FloatToStr(Laplas.egddata.myflmod[0].roughness);
         case Laplas.egddata.myflmod[0].ipowerroughness of
           1 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=0;
               end;
           2 : begin
                  FormSmagorinsky.CBpowermodel.ItemIndex:=1;
               end;
         end;
      end
      else
      begin
         FormSmagorinsky.Panrougness.Visible:=false;
      end;
      // поправка для течений с кривизной линий тока.
      FormSmagorinsky.CBSwirlAmendment.Checked:=Laplas.egddata.myflmod[0].bSwirlamendment;
      if (FormSmagorinsky.CBSwirlAmendment.Checked) then
      begin
         FormSmagorinsky.PanRichardson.Visible:=true;
         FormSmagorinsky.ERimult.Text:=FloatToStr(Laplas.egddata.myflmod[0].rRimult);
      end
        else
      begin
         FormSmagorinsky.PanRichardson.Visible:=false;
      end;
      // Установка модели Selective Smagorinsky
      FormSmagorinsky.CBSelectiveSmagorinsky.Checked:=Laplas.egddata.myflmod[0].bSelectiveSmagorinsky;
      if (FormSmagorinsky.CBSelectiveSmagorinsky.Checked) then
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=true;
         FormSmagorinsky.RGtypefiltr.ItemIndex:=Laplas.egddata.myflmod[0].itypefiltr;
         FormSmagorinsky.Eangle.Text:=FloatToStr(Laplas.egddata.myflmod[0].rSelectiveAngle);
      end
       else
      begin
         FormSmagorinsky.PSelectiveSmag.Visible:=false;
      end;
      FormSmagorinsky.ShowModal;
   end;
end;

end.
