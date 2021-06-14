unit Unitamgmanager;
// 26 октября 2016.
// Данный модуль предназначен для преодоления
// расходимости в проблемных задачах.
// Наиболее сильные настройки расположены в верхней части.
// В конце расположены настройки обращение к которым редко или
// которые имеют слабый эффект.
// Эти настройки предназначены только для домашнего кода РУМБА v.0.14
// когда он испытывает проблемы со сходимостью.
// Когда я пойму универсальные настройки на все случаи жизни то может быть
// я откажусь от этого меню. Но чтобы понять закономерности надо экспериментировать
// во всём диапазоне задач от самых слабых до самых жёстких выбирая те или иные
// настройки.
// Настройки вынесены в меню ещё и потому что компиляция после подсоединения массы
// библиотек занимает очень много времени и поэтому для быстрого экспериментирования
// сделано данное меню.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls, Vcl.AppEvnts;

type
  TForm_amg_manager = class(TForm)
    Panel1: TPanel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Buttondefault: TButton;
    Label7: TLabel;
    Label8: TLabel;
    Labelmagic: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    Label13: TLabel;
    Label14: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    ApplicationEvents1: TApplicationEvents;
    Label17: TLabel;
    Label18: TLabel;
    Label19: TLabel;
    Label20: TLabel;
    Label21: TLabel;
    CheckBoxprintlogTemperature: TCheckBox;
    CheckBoxprintlogSpeed: TCheckBox;
    CheckBoxprintlogPressure: TCheckBox;
    Label26: TLabel;
    LabelStress: TLabel;
    CheckBoxprintlogStress: TCheckBox;
    GroupBox1: TGroupBox;
    ComboBoxCFalgorithmandDataStruct_Temperature: TComboBox;
    ComboBoxCFalgorithmandDataStruct_Speed: TComboBox;
    ComboBoxCFalgorithmandDataStruct_Pressure: TComboBox;
    ComboBoxCFalgorithmandDataStruct_Stress: TComboBox;
    Label22: TLabel;
    CheckBoxTemperatureMatrixPortrait: TCheckBox;
    CheckBoxSpeedMatrixPortrait: TCheckBox;
    CheckBoxPressureMatrixPortrait: TCheckBox;
    CheckBoxStressMatrixPortrait: TCheckBox;
    ComboBoxSort: TComboBox;
    CheckBoxDiagonalDominant: TCheckBox;
    CheckBoxStrongTranspose: TCheckBox;
    PanelTemperature1: TPanel;
    Editthreshold: TEdit;
    ComboBoxmaximumreducedlevels: TComboBox;
    ComboBoxnFinnest: TComboBox;
    ComboBoxnumberpresmothers: TComboBox;
    ComboBoxnumberpostsweeps: TComboBox;
    ComboBoxmemorysize: TComboBox;
    ComboBoxinterpolation: TComboBox;
    PanelSpeed1: TPanel;
    EditthresholdSpeed: TEdit;
    ComboBoxmaximumreducedlevelsSpeed: TComboBox;
    ComboBoxnFinnestSpeed: TComboBox;
    ComboBoxnumberpresmothersSpeed: TComboBox;
    ComboBoxnumberpostsweepsSpeed: TComboBox;
    ComboBoxmemorysizeSpeed: TComboBox;
    ComboBoxInterpolationSpeed: TComboBox;
    PanelPressure1: TPanel;
    EditthresholdPressure: TEdit;
    ComboBoxmaximumreducedlevelsPressure: TComboBox;
    ComboBoxnFinnestPressure: TComboBox;
    ComboBoxnumberpresmothersPressure: TComboBox;
    ComboBoxnumberpostsweepsPressure: TComboBox;
    ComboBoxmemorysizePressure: TComboBox;
    ComboBoxinterpolationPressure: TComboBox;
    PanelStress1: TPanel;
    EditthresholdStress: TEdit;
    ComboBoxmaximumreducedlevelsStress: TComboBox;
    ComboBoxnFinnestStress: TComboBox;
    ComboBoxnumberpresmoothersStress: TComboBox;
    ComboBoxnumberpostsweepsStress: TComboBox;
    ComboBoxmemorysizeStress: TComboBox;
    ComboBoxinterpollationStress: TComboBox;
    PanelTemperature2: TPanel;
    CheckBoxtruncationT: TCheckBox;
    Edit_truncation_T: TEdit;
    EditmagicT: TEdit;
    ComboBoxsmoothertypeTemperature: TComboBox;
    ComboBoxcoarseningTemp: TComboBox;
    ComboBoxStabilizationTemp: TComboBox;
    PanelSpeed2: TPanel;
    CheckBoxtruncationSpeed: TCheckBox;
    Edit_truncation_Speed: TEdit;
    EditmagicSpeed: TEdit;
    ComboBoxsmoothertypeSpeed: TComboBox;
    ComboBoxcoarseningSpeed: TComboBox;
    ComboBoxStabilizationSpeed: TComboBox;
    PanelPressure2: TPanel;
    CheckBoxtruncationPressure: TCheckBox;
    Edit_truncation_Pressure: TEdit;
    EditmagicPressure: TEdit;
    ComboBoxsmoothertypePressure: TComboBox;
    ComboBoxcoarseningPressure: TComboBox;
    ComboBoxStabilizationPressure: TComboBox;
    PanelStress2: TPanel;
    CheckBoxtruncationStress: TCheckBox;
    Edittruncation_Stress: TEdit;
    EditmagicStress: TEdit;
    ComboBoxsmoothertypeStress: TComboBox;
    ComboBoxcoarseningStress: TComboBox;
    ComboBoxstabilizationStress: TComboBox;
    procedure ButtondefaultClick(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure CheckBoxtruncationTClick(Sender: TObject);
    procedure CheckBoxtruncationSpeedClick(Sender: TObject);
    procedure CheckBoxtruncationPressureClick(Sender: TObject);
    procedure Edittruncation_StressClick(Sender: TObject);
    procedure CheckBoxtruncationStressClick(Sender: TObject);
    procedure ComboBoxcoarseningTempChange(Sender: TObject);
    procedure ComboBoxcoarseningSpeedChange(Sender: TObject);
    procedure ComboBoxcoarseningPressureChange(Sender: TObject);
    procedure ComboBoxcoarseningStressChange(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form_amg_manager: TForm_amg_manager;

implementation

{$R *.dfm}

// запрет форме сворачиваться.
procedure TForm_amg_manager.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
     if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

 // Возвращает настройки решателя по умолчанию.
procedure TForm_amg_manager.ButtondefaultClick(Sender: TObject);
begin

     // 1 - QUICK SORT
     // 2 - HEAP SORT
     // 3 - TIM SORT (Тим Петерсон)
     // 4 - in development sort
     ComboBoxSort.ItemIndex:=0; // 0 - COUNTING SORT

     CheckBoxTemperatureMatrixPortrait.Checked:=false;
     CheckBoxSpeedMatrixPortrait.Checked:=false;
     CheckBoxPressureMatrixPortrait.Checked:=false;
     CheckBoxStressMatrixPortrait.Checked:=false;

    // 0 - AVL Tree, 1 - SPLAY Tree, 2 - Binary Heap.
    // 3 - Treap.
    ComboBoxCFalgorithmandDataStruct_Temperature.ItemIndex:=2;
    ComboBoxCFalgorithmandDataStruct_Speed.ItemIndex:=2;
    ComboBoxCFalgorithmandDataStruct_Pressure.ItemIndex:=2;
    ComboBoxCFalgorithmandDataStruct_Stress.ItemIndex:=2;

    // печать лога:
    CheckBoxprintlogTemperature.Checked:=true;
    CheckBoxprintlogSpeed.Checked:=true;
    CheckBoxprintlogPressure.Checked:=true;
    CheckBoxprintlogStress.Checked:=true;

    CheckBoxStrongTranspose.Checked:=true;

    // При возникновении проблем со сходимостью надо пробовать менять настройки
    // по умолчанию. Наиболее сильные настройки перечислены сверху.
    // Далее всё ниже и ниже эфеект изменения настроек падает.


    // Возвращает настройки по умолчанию.
    ComboBoxmaximumreducedlevels.ItemIndex:=0; // начиная с этого уровня вышестоящие уровни редуцируются (игнорируются).
    ComboBoxmaximumreducedlevelsSpeed.ItemIndex:=0;
    ComboBoxmaximumreducedlevelsPressure.ItemIndex:=0;
    ComboBoxmaximumreducedlevelsStress.ItemIndex:=0;

    ComboBoxinterpolation.ItemIndex:=3; // номер процедуры интерполляции.
    ComboBoxinterpolationSpeed.ItemIndex:=3;
    ComboBoxinterpolationPressure.ItemIndex:=3;
    ComboBoxinterpollationStress.ItemIndex:=3;

    ComboBoxnFinnest.ItemIndex:=1; // nFinnest
    ComboBoxnFinnestSpeed.ItemIndex:=1;
    ComboBoxnFinnestPressure.ItemIndex:=1;
    ComboBoxnFinnestStress.ItemIndex:=1;

    ComboBoxnumberpresmothers.ItemIndex:=1; // nu1
    ComboBoxnumberpresmothersSpeed.ItemIndex:=1;
    ComboBoxnumberpresmothersPressure.ItemIndex:=1;
    ComboBoxnumberpresmoothersStress.ItemIndex:=1;

    ComboBoxnumberpostsweeps.ItemIndex:=2; //nu2
    ComboBoxnumberpostsweepsSpeed.ItemIndex:=2;
    ComboBoxnumberpostsweepsPressure.ItemIndex:=2;
    ComboBoxnumberpostsweepsStress.ItemIndex:=2;

    ComboBoxmemorysize.ItemIndex:=5;  // memory number A size
    ComboBoxmemorysizeSpeed.ItemIndex:=5;
    ComboBoxmemorysizePressure.ItemIndex:=5;
    ComboBoxmemorysizeStress.ItemIndex:=5;


   // По умолчанию мы не используем ILU2 сглаживатель.
   // 0 - standart проекционный метод, 1 - ILU0, 1 - ILU2.
   ComboBoxsmoothertypeTemperature.ItemIndex:=0; // 0 - Standart.
   ComboBoxsmoothertypeSpeed.ItemIndex:=0; // 0 - Standart.
   ComboBoxsmoothertypePressure.ItemIndex:=0; // 0 - Standart.
   ComboBoxsmoothertypeStress.ItemIndex:=0; // 0 - Standart.

   if (FormatSettings.DecimalSeparator=',') then
   begin
      Editthreshold.Text:='0,24';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       Editthreshold.Text:='0.24';
   end;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      EditmagicT.Text:='0,4';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditmagicT.Text:='0.4';
   end;

   if (FormatSettings.DecimalSeparator=',') then
   begin
      EditmagicSpeed.Text:='0,4';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditmagicSpeed.Text:='0.4';
   end;

   if (FormatSettings.DecimalSeparator=',') then
   begin
      EditmagicPressure.Text:='0,4';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditmagicPressure.Text:='0.4';
   end;

   if (FormatSettings.DecimalSeparator=',') then
   begin
      EditmagicStress.Text:='0,4';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditmagicStress.Text:='0.4';
   end;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      EditthresholdSpeed.Text:='0,24';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditthresholdSpeed.Text:='0.24';
   end;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      EditthresholdPressure.Text:='0,24';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditthresholdPressure.Text:='0.24';
   end;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      EditthresholdStress.Text:='0,24';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       EditthresholdStress.Text:='0.24';
   end;

   // AmgSplitting standart vs RS2. 2.01.2017
   // ST option - 8.01.2017.
   // neg con only begin 19.01.2017
   ComboBoxcoarseningTemp.ItemIndex:=2; // standart coarsening  all con
   ComboBoxcoarseningSpeed.ItemIndex:=2; // standart coarsening all con
   ComboBoxcoarseningPressure.ItemIndex:=2; // standart coarsening all con
   ComboBoxcoarseningStress.ItemIndex:=2; // standart coarsening all con

   // Stabilization BiCGStab + amg (РУМБА).
   // 8.01.2017
   ComboBoxStabilizationTemp.ItemIndex:=0; // none
   ComboBoxStabilizationSpeed.ItemIndex:=0; // none
   ComboBoxStabilizationPressure.ItemIndex:=0; // none
   ComboBoxStabilizationStress.ItemIndex:=0; // none

    if (FormatSettings.DecimalSeparator=',') then
   begin
      Edit_truncation_T.Text:='0,2';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       Edit_truncation_T.Text:='0.2';
   end;
   Edit_truncation_T.Visible:=false;
   CheckBoxtruncationT.Checked:=true;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      Edit_truncation_Speed.Text:='0,2';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       Edit_truncation_Speed.Text:='0.2';
   end;
   Edit_truncation_Speed.Visible:=false;
   CheckBoxtruncationSpeed.Checked:=true;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      Edit_truncation_Pressure.Text:='0,2';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       Edit_truncation_Pressure.Text:='0.2';
   end;
   Edit_truncation_Pressure.Visible:=false;
   CheckBoxtruncationPressure.Checked:=true;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      Edittruncation_Stress.Text:='0,2';
   end
    else  if (FormatSettings.DecimalSeparator='.') then
   begin
       Edittruncation_Stress.Text:='0.2';
   end;
   Edittruncation_Stress.Visible:=false;
   CheckBoxtruncationStress.Checked:=true;

end;

// truncation interpolation for Pressure click.
procedure TForm_amg_manager.CheckBoxtruncationPressureClick(Sender: TObject);
var
  s1 : String;
  i : Integer;
begin
   if (CheckBoxtruncationPressure.Checked) then
   begin
      Edit_truncation_Pressure.Visible:=false;
   end
   else
   begin
      Edit_truncation_Pressure.Visible:=true;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=Edit_truncation_Pressure.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]=',') then s1[i]:='.';
         end;
         Edit_truncation_Pressure.Text:=s1;
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=Edit_truncation_Pressure.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]='.') then s1[i]:=',';
         end;
         Edit_truncation_Pressure.Text:=s1;
      end;
   end;
end;

// truncation interpolation for Speed click.
procedure TForm_amg_manager.CheckBoxtruncationSpeedClick(Sender: TObject);
var
  s1 : String;
  i : Integer;
begin
   if (CheckBoxtruncationSpeed.Checked) then
   begin
      Edit_truncation_Speed.Visible:=false;
   end
   else
   begin
      Edit_truncation_Speed.Visible:=true;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=Edit_truncation_Speed.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]=',') then s1[i]:='.';
         end;
         Edit_truncation_Speed.Text:=s1;
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=Edit_truncation_Speed.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]='.') then s1[i]:=',';
         end;
         Edit_truncation_Speed.Text:=s1;
      end;
   end;
end;

procedure TForm_amg_manager.CheckBoxtruncationStressClick(Sender: TObject);
var
   s1 : String;
   i : Integer;
begin
   if (CheckBoxtruncationStress.Checked) then
   begin
      Edittruncation_Stress.Visible:=false;
   end
   else
   begin
      Edittruncation_Stress.Visible:=true;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=Edittruncation_Stress.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]=',') then s1[i]:='.';
         end;
         Edittruncation_Stress.Text:=s1;
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=Edittruncation_Stress.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]='.') then s1[i]:=',';
         end;
         Edittruncation_Stress.Text:=s1;
      end;
   end;
end;

// truncation interpolation for Temperature click.
procedure TForm_amg_manager.CheckBoxtruncationTClick(Sender: TObject);
var
  s1 : String;
  i : Integer;
begin
   if (CheckBoxtruncationT.Checked) then
   begin
      Edit_truncation_T.Visible:=false;
   end
   else
   begin
      Edit_truncation_T.Visible:=true;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=Edit_truncation_T.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]=',') then s1[i]:='.';
         end;
         Edit_truncation_T.Text:=s1;
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=Edit_truncation_T.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]='.') then s1[i]:=',';
         end;
         Edit_truncation_T.Text:=s1;
      end;
   end;
end;




procedure TForm_amg_manager.ComboBoxcoarseningPressureChange(Sender: TObject);
begin
    if ((ComboBoxcoarseningPressure.ItemIndex=8)or(ComboBoxcoarseningPressure.ItemIndex=9)) then
   begin
      // PMIS, HMIS
      CheckBoxtruncationPressure.Checked:=false; // Применяем усечение интерполяции при PMIS или HMIS.
      // threshold зависит от номера уровня и эта зависимость задана программистом.
      EditthresholdPressure.Visible:=false; // creator depend
      if (ComboBoxcoarseningPressure.ItemIndex=8) then
      begin
         // PMIS
         // Используется собственная специальная процедура интерполяции определенная программистом.
         //ComboBoxinterpolationPressure.Visible:=true;//false; // creator depend
      end
      else
      begin
        // ComboBoxinterpolationPressure.Visible:=true;
      end;
      // в HMIS на первом уровне выбор процедуры интекрполяции всё же используется.
   end
   else
   begin
      EditthresholdPressure.Visible:=true;
     // ComboBoxinterpolationPressure.Visible:=true;
   end;
end;

procedure TForm_amg_manager.ComboBoxcoarseningSpeedChange(Sender: TObject);
begin
    if ((ComboBoxcoarseningSpeed.ItemIndex=8)or(ComboBoxcoarseningSpeed.ItemIndex=9)) then
   begin
      // PMIS, HMIS
      CheckBoxtruncationSpeed.Checked:=false; // Применяем усечение интерполяции при PMIS или HMIS.
      // threshold зависит от номера уровня и эта зависимость задана программистом.
      EditthresholdSpeed.Visible:=false; // creator depend
      if (ComboBoxcoarseningSpeed.ItemIndex=8) then
      begin
         // PMIS
         // Используется собственная специальная процедура интерполяции определенная программистом.
        // ComboBoxinterpolationSpeed.Visible:=true;//false; // creator depend
      end
      else
      begin
         //ComboBoxinterpolationSpeed.Visible:=true;
      end;
      // в HMIS на первом уровне выбор процедуры интекрполяции всё же используется.
   end
   else
   begin
      EditthresholdSpeed.Visible:=true;
      //ComboBoxinterpolationSpeed.Visible:=true;
   end;
end;

procedure TForm_amg_manager.ComboBoxcoarseningStressChange(Sender: TObject);
begin
   if ((ComboBoxcoarseningStress.ItemIndex=8)or(ComboBoxcoarseningStress.ItemIndex=9)) then
   begin
      // PMIS, HMIS
      CheckBoxtruncationStress.Checked:=false; // Применяем усечение интерполяции при PMIS или HMIS.
      // threshold зависит от номера уровня и эта зависимость задана программистом.
      EditthresholdStress.Visible:=false; // creator depend
      if (ComboBoxcoarseningStress.ItemIndex=8) then
      begin
         // PMIS
         // Используется собственная специальная процедура интерполяции определенная программистом.
       //  ComboBoxinterpollationStress.Visible:=true;//false; // creator depend
      end
      else
      begin
         // HMIS
        // ComboBoxinterpollationStress.Visible:=true;
      end;
      // в HMIS на первом уровне выбор процедуры интекрполяции всё же используется.
   end
   else
   begin
      EditthresholdStress.Visible:=true;
     // ComboBoxinterpollationStress.Visible:=true;
   end;
end;

procedure TForm_amg_manager.ComboBoxcoarseningTempChange(Sender: TObject);
begin
   // Temperature
   if ((ComboBoxcoarseningTemp.ItemIndex=8)or(ComboBoxcoarseningTemp.ItemIndex=9)) then
   begin
      // PMIS, HMIS
      ComboBoxmemorysize.ItemIndex:=0;// size ierarhion 4A
      CheckBoxtruncationT.Checked:=false; // Применяем усечение интерполяции при PMIS или HMIS.
      // threshold зависит от номера уровня и эта зависимость задана программистом.
      Editthreshold.Visible:=false; // creator depend
      if (ComboBoxcoarseningTemp.ItemIndex=8) then
      begin
         // PMIS
         // Используется собственная специальная процедура интерполяции определенная программистом.
        // ComboBoxinterpolation.Visible:=true;//false; // creator depend
      end
      else
      begin
         // HMIS
        // ComboBoxinterpolation.Visible:=true;
      end;
      // в HMIS на первом уровне выбор процедуры интекрполяции всё же используется.
   end
   else
   begin
      Editthreshold.Visible:=true;
     // ComboBoxinterpolation.Visible:=true;
   end;
end;

procedure TForm_amg_manager.Edittruncation_StressClick(Sender: TObject);
var
  s1 : String;
  i : Integer;
begin
    if (CheckBoxtruncationStress.Checked) then
   begin
      Edittruncation_Stress.Visible:=false;
   end
   else
   begin
      Edittruncation_Stress.Visible:=true;
      if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=Edittruncation_Stress.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]=',') then s1[i]:='.';
         end;
         Edittruncation_Stress.Text:=s1;
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=Edittruncation_Stress.Text;
         for i := 1 to length(s1) do
         begin
            if (s1[i]='.') then s1[i]:=',';
         end;
         Edittruncation_Stress.Text:=s1;
      end;
   end;
end;

end.
