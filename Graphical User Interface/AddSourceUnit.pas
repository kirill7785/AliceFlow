unit AddSourceUnit;
// изменение свойств источника тепла

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, ExtCtrls, StdCtrls;

type
  TAddSourceForm = class(TForm)
    Panelglobalcontainer: TPanel;
    RadioGroup1: TRadioGroup;
    Bapply: TButton;
    Panelinfo: TPanel;
    Lname: TLabel;
    Ename: TEdit;
    PanelGeometry: TPanel;
    GroupBox1: TGroupBox;
    LxS: TLabel;
    LyS: TLabel;
    LzS: TLabel;
    LxE: TLabel;
    LyE: TLabel;
    LzE: TLabel;
    ExS: TEdit;
    EyS: TEdit;
    EzS: TEdit;
    ExE: TEdit;
    EyE: TEdit;
    EzE: TEdit;
    RadioGroupPlane: TRadioGroup;
    PanelProperties: TPanel;
    GBpowerdef: TGroupBox;
    Label1: TLabel;
    LW: TLabel;
    RGpowertype: TRadioGroup;
    Ptempdefloc: TPanel;
    Label2: TLabel;
    Label3: TLabel;
    CBtableid: TComboBox;
    EOperoffsetdrain: TEdit;
    Epower: TEdit;
    procedure BapplyClick(Sender: TObject);
    procedure RadioGroupPlaneClick(Sender: TObject);
    procedure RGpowertypeClick(Sender: TObject);
    procedure RadioGroup1Click(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  AddSourceForm: TAddSourceForm;

implementation
uses
     VisualUnit, UnitVariables;
{$R *.dfm}

procedure TAddSourceForm.BapplyClick(Sender: TObject);
var
   k : Integer;
   // вспомогательные переменные для обработки исключительной ситуации
   bOk : Boolean;
   s1, s2, s3, s4, s5, s6, spow : String;
   r1, r2, r3, r4, r5, r6, rpow : Real;
   buf : TPlane;

begin
   // инициализация :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;

   s1:=Trim(Epower.Text);
   (*
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
   *)
   if (FormatSettings.DecimalSeparator=',') then
      begin
         s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
      end;
      if (FormatSettings.DecimalSeparator='.') then
      begin
         s1:=StringReplace(s1,',','.',[rfReplaceAll]);
      end;
   Epower.Text:=s1;


   s1:=Trim(ExS.Text);
   (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
       s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
    end;
   ExS.Text:=s1;

    s1:=Trim(EyS.Text);
    (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
      s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
    end;
   EyS.Text:=s1;

    s1:=Trim(EzS.Text);
    (*
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
   *)
    if (FormatSettings.DecimalSeparator='.') then
    begin
       s1:=StringReplace(s1,',','.',[rfReplaceAll]);
    end;

    if (FormatSettings.DecimalSeparator=',') then
    begin
       s1:=StringReplace(s1,'.',',',[rfReplaceAll]);
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

    s1:=Trim(EOperoffsetdrain.Text);
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
   EOperoffsetdrain.Text:=s1;


   // ввод данных об источнике
   k:=Laplas.itek;
   //with (Laplas.source[k]) do
   buf:=Laplas.sourcepublic[k];
  //with (buf) do
   //begin
      bOk:=true; // признак правельности ввода
      buf.iPlane:=RadioGroupPlane.ItemIndex+1;
      spow:=Epower.Text;
      if bOk then rpow:=FormVariables.my_real_convert(spow,bOk);
      // тип задания мощности : 0 - константа, 1 - таблично.
      buf.itempdep:=RGpowertype.ItemIndex;
      if (buf.itempdep=1) then
      begin
         buf.soperatingoffsetdrain:=EOperoffsetdrain.Text;
         if bOk then buf.operatingoffsetdrain:=FormVariables.my_real_convert(buf.soperatingoffsetdrain,bOk);
         buf.id_table:=CBtableid.ItemIndex; // уникальный номер таблицы
      end;

      buf.name:=Ename.Text; // имя элемента
      // корректировка имени объекта чтобы избежать совпадающих имён.
      Laplas.correctobjname('s',buf.name,k);
      Ename.Text:=buf.name;


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
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учётом подстановки
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
               if bOk then r3:=FormVariables.my_real_convert(s3,bOk);  // с учётом подстановки
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
               if bOk then r3:=r1; // с учётом подстановки
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
         buf.spower:=spow;

         buf.xS:=r1;  // числовые размеры
         buf.yS:=r2;  // заданные пользователем
         buf.xE:=r3;  // с учётом подстановки
         buf.yE:=r4;  // значений переменных.
         buf.zS:=r5;
         buf.zE:=r6;
         buf.Power:=rpow;
      end;

   //end;

   Laplas.sourcepublic[k]:=buf;
   with Laplas do
   begin
      ReadyPaint;
   end;

end;

procedure TAddSourceForm.RadioGroup1Click(Sender: TObject);
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
            PanelInfo.Visible:=false;
            PanelGeometry.Visible:=true;
            PanelProperties.Visible:=false;
         end;
     2 : begin
            // Properties
            PanelInfo.Visible:=false;
            PanelGeometry.Visible:=false;
            PanelProperties.Visible:=true;
         end;
   end;
end;

procedure TAddSourceForm.RadioGroupPlaneClick(Sender: TObject);
begin
   // смена плоскости в которой лежит источник тепла
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

procedure TAddSourceForm.RGpowertypeClick(Sender: TObject);
var
    i : Integer;
begin
    // смена способа задания мощности :
    // мощность может задаваться либо постоянной (константа)
    // либо в виде таблицы как зависимость от температуры и смещения стока.
    case RGpowertype.ItemIndex of
       0 : // const
           begin
              Ptempdefloc.Visible:=false;
              Label1.Caption:='power';
              LW.Caption:='W';
           end;
       1 : // power define
           begin
              if (Laplas.iltdp=0) then
              begin
                  Application.MessageBox('Plese Define -> Power Table create','Please define',MB_OK);
                  Laplas.MainMemo.Lines.Add('Plese Define -> Power Table create');
                  RGpowertype.ItemIndex:=0; // постоянная мощность.
                  Ptempdefloc.Visible:=false;
                  Label1.Caption:='power';
                  LW.Caption:='W';
              end
               else
              begin
                 // обновление списка доступных таблиц.
                 CBtableid.Clear;
                 for i:=0 to Laplas.iltdp-1 do
                 begin
                    CBtableid.AddItem(IntToStr(i),Sender);
                 end;
                 CBtableid.ItemIndex:=Laplas.source[Laplas.itek].id_table;  // уникальный номер таблично заданной мощности.
                 EOperoffsetdrain.Text:=Laplas.source[Laplas.itek].soperatingoffsetdrain; // напряжение на стоке
                 Ptempdefloc.Visible:=true;
                 // mult power - это множитель
                 // на который домножается реальная мощность.
                 // Это удобно если считается только половина реального источника
                 // в симметричной постановке, тогда мощность надо разделить пополам
                 // что соответствует mult power равным 0.5.
                 Label1.Caption:='mult power';
                 LW.Caption:='';
              end;
           end;
    end;
end;

end.

