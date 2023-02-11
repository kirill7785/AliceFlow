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
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants,
  System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ExtCtrls, Vcl.StdCtrls;

type
  TFormRadiation = class(TForm)
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
    ButtonApply: TButton;
    CheckBoxinternalRadiation: TCheckBox;
    procedure FormCreate(Sender: TObject);
    procedure ButtonApplyClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }

  end;

var
  FormRadiation: TFormRadiation;

implementation

{$R *.dfm}
uses
    Math, VisualUnit,  Unitrectangularplot, UnitVariables;


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
       // Закрываем экранную форму после успешного ввода.
       Close;
    end
     else
    begin
       ShowMessage('Your input is incorrect. Please, reinput...');
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
end;

end.

