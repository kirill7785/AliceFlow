unit AddWallUnit;
// в этом модуле добавл€етс€ элемент стенка
// стенка может лежать только на границе кабинета.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TAddWallForm = class(TForm)
    Panelglobalcontainer: TPanel;
    RadioGroup1: TRadioGroup;
    PanelInfo: TPanel;
    Ename: TEdit;
    Lname: TLabel;
    Bapply: TButton;
    PanelGeometry: TPanel;
    RadioGroupPlane: TRadioGroup;
    GroupBoxSize: TGroupBox;
    LxS: TLabel;
    LyS: TLabel;
    LzS: TLabel;
    LxE: TLabel;
    LyE: TLabel;
    LzE: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    ExS: TEdit;
    EyS: TEdit;
    EzS: TEdit;
    ExE: TEdit;
    EyE: TEdit;
    EzE: TEdit;
    PanelProperties: TPanel;
    GroupBoxFLOW: TGroupBox;
    RadioGroupflowtype: TRadioGroup;
    GroupBoxvelcomp: TGroupBox;
    LabelVx: TLabel;
    LabelVy: TLabel;
    LabelVz: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    EditVx: TEdit;
    EditVy: TEdit;
    EditVz: TEdit;
    GroupBoxpressure: TGroupBox;
    Labelpressure: TLabel;
    Label1: TLabel;
    Editpress: TEdit;
    GroupBoxtemper: TGroupBox;
    RadioGroupBonConTemp: TRadioGroup;
    PaneltemperatureBC: TPanel;
    LTemp: TLabel;
    Label5: TLabel;
    Etemp: TEdit;
    Panelemissivity: TPanel;
    Label12: TLabel;
    Editemissivity: TEdit;
    GroupBoxThermalStress: TGroupBox;
    Labeldeformation: TLabel;
    ComboBoxDeformationBoundaryConditon: TComboBox;
    EditForce: TEdit;
    LabelForce: TLabel;
    procedure BapplyClick(Sender: TObject);
    procedure RadioGroupPlaneClick(Sender: TObject);
    procedure RadioGroupflowtypeClick(Sender: TObject);
    procedure RadioGroupBonConTempClick(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
    procedure ComboBoxDeformationBoundaryConditonChange(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  AddWallForm: TAddWallForm;

implementation

uses VisualUnit, UnitVariables;

{$R *.dfm}

// редактирование свойств стенки
procedure TAddWallForm.BapplyClick(Sender: TObject);
var
    k : Integer;
    // вспомогательные переменные дл€ обработки исключительной ситуации
    bOk : Boolean;
    s1, s2, s3, s4, s5, s6 : String;
    r1, r2, r3, r4, r5, r6 : Real;
    buf : TPlane;

begin
   // инициализаци€ :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;

   s1:=Trim(Etemp.Text);
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
   Etemp.Text:=s1;


   s1:=Trim(EditVx.Text);
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
   EditVx.Text:=s1;

   s1:=Trim(EditVy.Text);
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
   EditVy.Text:=s1;

   s1:=Trim(EditVz.Text);
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
   EditVz.Text:=s1;

   s1:=Trim(ExS.Text);
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
   ExS.Text:=s1;

    s1:=Trim(EyS.Text);
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
   EyS.Text:=s1;

   s1:=Trim(EzS.Text);
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
   EzS.Text:=s1;

   s1:=Trim(ExE.Text);
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
   ExE.Text:=s1;

    s1:=Trim(EyE.Text);
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
   EyE.Text:=s1;

   s1:=Trim(EzE.Text);
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
   EzE.Text:=s1;


    // ввод данных о стенке
   k:=Laplas.itek;
   buf:=Laplas.wallpublic[k];

   //with (buf) do
   //begin
      buf.iPlane:=RadioGroupPlane.ItemIndex+1; // плоскость в которой лежит стенка
      buf.family:=RadioGroupBonConTemp.ItemIndex+1; // род граничных условий дл€ температуры.
      if ((buf.family=1)or(buf.family=3)or(buf.family=4)) then buf.Tamb:=StrToFloat(Etemp.Text)  // посто€нна€  температура
      else buf.Tamb:=0.0;

       s1:=Trim(Editemissivity.Text);
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
       Editemissivity.Text:=s1;
       bOk:=true;
      if (buf.family=4) then
      begin
         //buf.emissivity:=StrToFloat(Editemissivity.Text);
         buf.semissivity:=Trim(Editemissivity.Text);
         buf.emissivity:=FormVariables.my_real_convert(Editemissivity.Text,bOk);
      end;
      if (buf.family=3) then
      begin
         //buf.heat_transfer_coefficient:=StrToFloat(Editemissivity.Text);
         buf.sheat_transfer_coefficient:=Trim(Editemissivity.Text);
         buf.heat_transfer_coefficient:=FormVariables.my_real_convert(Editemissivity.Text,bOk);
      end;
      buf.HF:=0.0; // только нулевой тепловой поток
      buf.name:=Trim(Ename.Text); // им€ элемента
      // корректировка имени объекта чтобы избежать совпадающих имЄн.
      Laplas.correctobjname('w',buf.name,Laplas.itek);
      Ename.Text:=buf.name;

      bOk:=true;

      buf.ithermal_stress_boundary_condition:= ComboBoxDeformationBoundaryConditon.ItemIndex;
      s1:=Trim(EditForce.Text);
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
       EditForce.Text:=s1;
       if (buf.ithermal_stress_boundary_condition=8) then
          begin
             buf.xForce:=StrToFloat(EditForce.Text);
          end;
          if (buf.ithermal_stress_boundary_condition=9) then
          begin
             buf.yForce:=StrToFloat(EditForce.Text);
          end;
          if (buf.ithermal_stress_boundary_condition=10) then
          begin
             buf.zForce:=StrToFloat(EditForce.Text);
          end;

      case buf.iPlane of
        1 : // XY
            begin
               s1:=ExS.Text;  // параметризованные
               s2:=EyS.Text;  // геометрические
               s3:=ExE.Text;  // размеры
               s4:=EyE.Text;  // заданные
               s5:=EzS.Text;  // пользователем
               s6:=s5;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учЄтом подстановки
               if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // значений переменных.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=r5;
            end;
        2 : // XZ
            begin
               s1:=ExS.Text; // параметризованные
               s2:=EyS.Text; // геометрические
               s3:=ExE.Text; // размеры
               s4:=s2;      // заданные
               s5:=EzS.Text; // пользователем
               s6:=EzE.Text;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учЄтом подстановки
               if bOk then r4:=r2;  // значений переменных.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=FormVariables.my_real_convert(s6,bOk);
            end;
        3 : // YZ
            begin
               s1:=ExS.Text; // параметризованные
               s2:=EyS.Text; // геометрические
               s3:=s1;  // размеры
               s4:=EyE.Text; // заданные
               s5:=EzS.Text; // пользователем
               s6:=EzE.Text;

               if bOk then r1:=FormVariables.my_real_convert(s1,bOk);  // числовые размеры
               if bOk then r2:=FormVariables.my_real_convert(s2,bOk);  // заданные пользователем
               if bOk then r3:=r1; // с учЄтом подстановки
               if bOk then r4:=FormVariables.my_real_convert(s4,bOk);  // значений переменных.
               if bOk then r5:=FormVariables.my_real_convert(s5,bOk);
               if bOk then r6:=FormVariables.my_real_convert(s6,bOk);
            end;
      end; // case


      if (bOk) then
      begin
         buf.sxS:=s1;  // параметризованные
         buf.syS:=s2;  // геометрические
         buf.sxE:=s3;  // размеры
         buf.syE:=s4;  // заданные
         buf.szS:=s5;  // пользователем
         buf.szE:=s6;

         buf.xS:=r1;  // числовые размеры
         buf.yS:=r2;  // заданные пользователем
         buf.xE:=r3;  // с учЄтом подстановки
         buf.yE:=r4;  // значений переменных.
         buf.zS:=r5;
         buf.zE:=r6;
      end;

      // задание гидродинамических условий:
      case RadioGroupflowtype.ItemIndex of
        0 : begin
               //velocity
               buf.Vx:=StrToFloat(EditVx.Text);
               buf.Vy:=StrToFloat(EditVy.Text);
               buf.Vz:=StrToFloat(EditVz.Text);
               buf.P:=0.0;
               buf.bsymmetry:=false;
               buf.bpressure:=false;
               buf.bopening:=false;
            end;
        1 : begin
               // pressure
               buf.Vx:=0.0; buf.Vy:=0.0; buf.Vz:=0.0;
               buf.P:=StrToFloat(Editpress.Text);
               buf.bsymmetry:=false;
               buf.bpressure:=true;
               buf.bopening:=false;
            end;
        2 : begin
               // symmetry
               buf.Vx:=0.0; buf.Vy:=0.0; buf.Vz:=0.0;
               buf.P:=0.0;
               buf.bsymmetry:=true;
               buf.bpressure:=false;
               buf.bopening:=false;
            end;
        3 : begin
               // opening
               buf.Vx:=0.0; buf.Vy:=0.0; buf.Vz:=0.0;
               buf.P:=0.0;
               buf.bsymmetry:=false;
               buf.bpressure:=false;
               buf.bopening:=true;
            end;
       end;

   //end;
   Laplas.wallpublic[Laplas.itek]:=buf;
   //Laplas.wall[k]:=buf;
   with Laplas do
   begin
      ReadyPaint;
   end;
end;

procedure TAddWallForm.RadioGroupPlaneClick(Sender: TObject);
begin
   // смена плоскости в которой лежит твЄрда€ стенка
   case (RadioGroupPlane.ItemIndex+1) of
     1 : // XY
         begin
            ExE.Visible:=true;
            LxE.Visible:=true;
            EyE.Visible:=true;
            LyE.Visible:=true;
            EzE.Visible:=false;
            LzE.Visible:=false;
         end;
     2 : // XZ
         begin
            ExE.Visible:=true;
            LxE.Visible:=true;
            EyE.Visible:=false;
            LyE.Visible:=false;
            EzE.Visible:=true;
            LzE.Visible:=true;
         end;
     3 : // YZ
         begin
            ExE.Visible:=false;
            LxE.Visible:=false;
            EyE.Visible:=true;
            LyE.Visible:=true;
            EzE.Visible:=true;
            LzE.Visible:=true;
         end;
   end;
end;

procedure TAddWallForm.RadioGroupflowtypeClick(Sender: TObject);
var
    k : Integer;
begin
   k:=Laplas.itek;
   // реакци€ интерфейса на смену
   // граничного услови€ по скорости
   case RadioGroupflowtype.ItemIndex of
     0 : begin
            //velocity
            GroupBoxpressure.Visible:=false;
            GroupBoxvelcomp.Visible:=true;
         end;
     1 : begin
            // pressure
            GroupBoxpressure.Visible:=true;
            GroupBoxvelcomp.Visible:=false;
         end;
     2 : begin
            // symmetry
            GroupBoxpressure.Visible:=false;
            GroupBoxvelcomp.Visible:=false;
         end;
      3 : begin
            // opening
            GroupBoxpressure.Visible:=false;
            GroupBoxvelcomp.Visible:=false;
         end;
   end;

   with (Laplas.wall[k]) do
   begin
      EditVx.Text:=FloatToStr(Vx);
      EditVy.Text:=FloatToStr(Vy);
      EditVz.Text:=FloatToStr(Vz);
      Editpress.Text:=FloatToStr(P);
   end;
end;

procedure TAddWallForm.ComboBoxDeformationBoundaryConditonChange(
  Sender: TObject);
begin
      Laplas.wall[Laplas.itek].ithermal_stress_boundary_condition:= ComboBoxDeformationBoundaryConditon.ItemIndex;
      if (Laplas.wall[Laplas.itek].ithermal_stress_boundary_condition<8) then
      begin
          EditForce.Visible:=false;
          LabelForce.Visible:=false;
      end
       else
      begin
          EditForce.Visible:=true;
          LabelForce.Visible:=true;
          if (Laplas.wall[Laplas.itek].ithermal_stress_boundary_condition=8) then
          begin
             AddWallForm.EditForce.Text:=FloatToStr(Laplas.wall[Laplas.itek].xForce);
          end;
          if (Laplas.wall[Laplas.itek].ithermal_stress_boundary_condition=9) then
          begin
             AddWallForm.EditForce.Text:=FloatToStr(Laplas.wall[Laplas.itek].yForce);
          end;
          if (Laplas.wall[Laplas.itek].ithermal_stress_boundary_condition=10) then
          begin
             AddWallForm.EditForce.Text:=FloatToStr(Laplas.wall[Laplas.itek].zForce);
          end;
      end;
end;

procedure TAddWallForm.RadioGroup1Click(Sender: TObject);
begin
   case RadioGroup1.ItemIndex of
      0 : begin
         // Info
         PanelInfo.Visible:=true;
         PanelGeometry.Visible:=false;
         PanelProperties.Visible:=false;
      end;
      1 : begin
        // Geometry
        PanelGeometry.Visible:=true;
        PanelInfo.Visible:=false;
        PanelProperties.Visible:=false;
      end;
      2 : begin
        // Properties
        PanelGeometry.Visible:=false;
        PanelInfo.Visible:=false;
        PanelProperties.Visible:=true;
      end;
   end;
end;

procedure TAddWallForm.RadioGroupBonConTempClick(Sender: TObject);
var
    k : Integer;
begin
   k:=Laplas.itek;
   // смена типа граничного услови€ по температуре
   case  RadioGroupBonConTemp.ItemIndex of
     0 : begin
            // задана температура
            PaneltemperatureBC.Visible:=true;
            Etemp.Text:=FloatToStr(Laplas.wall[k].Tamb);
            Panelemissivity.Visible:=false;
         end;
     1 : begin
            // однородное условие Ќеймана
            PaneltemperatureBC.Visible:=false;
            Panelemissivity.Visible:=false;
         end;
     2 : begin
            // условие третьего рода Ќьютона-–ихмана
            Panelemissivity.Visible:=true;
            PaneltemperatureBC.Visible:=true;
            Etemp.Text:=FloatToStr(Laplas.wall[k].Tamb);
            Editemissivity.Text:=FloatToStr(Laplas.wall[k].heat_transfer_coefficient);
            Label12.Caption:='heat transfer coeff';
         end;
     3 : begin
            // ”словие —тефана-Ѕольцмана
            PaneltemperatureBC.Visible:=true;
            Panelemissivity.Visible:=true;
            Etemp.Text:=FloatToStr(Laplas.wall[k].Tamb);
            Editemissivity.Text:=FloatToStr(Laplas.wall[k].emissivity);
            Label12.Caption:='emissivity';
         end;
   end;
end;

end.
