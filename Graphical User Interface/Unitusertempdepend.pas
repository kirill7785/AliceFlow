unit Unitusertempdepend;
// Кусочно-линейное пользовательское задание параметров расчёта.
// Зависящие от температуры в градусах Цельсия.
// 17-18-19 ноября 2016, теплопроводность и теплоёмкость зависящие от температуры.
// 19-20 ноября 2016, объёмная мощность зависящая от температуры.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormusertempdepend = class(TForm)
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Memopiecewiseproperties: TMemo;
    Label5: TLabel;
    ButtonApply: TButton;
    ButtonView: TButton;
    Label6: TLabel;
    ComboBoxtemperatureUnit: TComboBox;
    procedure ButtonApplyClick(Sender: TObject);
    procedure ButtonViewClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    // 0 - cp,
    // 1 - lam,
    // 2 - power.
    identifier : Integer;
  end;

var
  Formusertempdepend: TFormusertempdepend;

implementation

{$R *.dfm}

uses VisualUnit, Unitrectangularplot, UnitVariables;

// Считывание табличной зависимости для теплопроводности или теплоёмкости
// материалов.
procedure TFormusertempdepend.ButtonApplyClick(Sender: TObject);
var
   s, sub : String;
   i,j : Integer;
   bOk : Boolean;
begin
   bOk:=true;
   j:=0;
   for i := 0 to Memopiecewiseproperties.Lines.Count-1 do
   begin
      s:= Memopiecewiseproperties.Lines[i];
      s:=Trim(s);
      if (length(s)>0) then
      begin
         sub:=Copy(s,1,Pos(' ',s)-1);
         sub:=Trim(sub);
         if (identifier=1) then
         begin
            // thermal conductivity
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam,j+1);
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam,j+1);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp,j+1);
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp,j+1);
         end;
         if (identifier=2) then
         begin
            // volume power
            SetLength(Laplas.body[Laplas.itek].temp_power,j+1);
            SetLength(Laplas.body[Laplas.itek].arr_power,j+1);
            SetLength(Laplas.body[Laplas.itek].arr_s_power,j+1);
         end;

         if (FormatSettings.DecimalSeparator='.') then
         begin
            // заменить все запятые в файле на точки.
            sub:=StringReplace(sub,',','.',[rfReplaceAll]);
         end;
         if (FormatSettings.DecimalSeparator=',') then
         begin
            // заменить все запятые в файле на точки.
            sub:=StringReplace(sub,'.',',',[rfReplaceAll]);
         end;
         if (identifier=1) then
         begin
            // thermal conductivity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]:=StrToFloat(sub);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]:=StrToFloat(sub);
         end;
         if (identifier=2) then
         begin
            // power
            Laplas.body[Laplas.itek].temp_power[j]:=StrToFloat(sub);
         end;

         s:=Copy(s,Pos(' ',s)+1,length(s));
         s:=Trim(s);
         if (FormatSettings.DecimalSeparator='.') then
         begin
            // заменить все запятые в файле на точки.
            s:=StringReplace(s,',','.',[rfReplaceAll]);
         end;
         if (FormatSettings.DecimalSeparator=',') then
         begin
            // заменить все запятые в файле на точки.
            s:=StringReplace(s,'.',',',[rfReplaceAll]);
         end;
         if (identifier=1) then
         begin
            // thermal conductivity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[j]:=StrToFloat(s);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[j]:=StrToFloat(s);
         end;
         if (identifier=2) then
         begin
            // power
            Laplas.body[Laplas.itek].arr_s_power[j]:=Trim(s);
            //Laplas.body[Laplas.itek].arr_power[j]:=StrToFloat(s);
            if (bOk) then
            begin
               Laplas.body[Laplas.itek].arr_power[j]:=FormVariables.my_real_convert(s,bOk);
            end;
         end;
         inc(j);
      end;
   end;
   if (identifier=1) then
   begin
      // thermal conductivity
      Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam:=j;
   end;
   if (identifier=0) then
   begin
      // heat capacity
      Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp:=j;
   end;
   if (identifier=2) then
   begin
      // power
      Laplas.body[Laplas.itek].n_power:=j;
   end;
   if (ComboBoxtemperatureUnit.ItemIndex=1) then
   begin
     // K 2 C
     if (identifier=1) then
     begin
        // thermal conductivity
        for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam-1 do
        begin
           Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]-273.0;
        end;
     end;
     if (identifier=0) then
     begin
        // heat capacity
        for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp-1 do
        begin
           Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]-273.0;
        end;
     end;
     if (identifier=2) then
     begin
        // power
        for j := 0 to Laplas.body[Laplas.itek].n_power-1 do
        begin
           Laplas.body[Laplas.itek].temp_power[j]:=Laplas.body[Laplas.itek].temp_power[j]-273.0;
        end;
     end;
   end;

end;

// Визуализация введенной пользователем зависимости.
procedure TFormusertempdepend.ButtonViewClick(Sender: TObject);
var
   s, sub : String;
   i,j : Integer;
   bOk : Boolean;
begin
   bOk:=true;
   // Этап 1 считывание пользовательской зависимости в память.
   j:=0;
   for i := 0 to Memopiecewiseproperties.Lines.Count-1 do
   begin
      s:= Memopiecewiseproperties.Lines[i];
      s:=Trim(s);
      if (length(s)>0) then
      begin
         sub:=Copy(s,1,Pos(' ',s)-1);
         sub:=Trim(sub);
         if (identifier=1) then
         begin
            // thermal conductivity
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam,j+1);
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam,j+1);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp,j+1);
            SetLength(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp,j+1);
         end;
         if (identifier=2) then
         begin
            // volume power
            SetLength(Laplas.body[Laplas.itek].temp_power,j+1);
            SetLength(Laplas.body[Laplas.itek].arr_power,j+1);
            SetLength(Laplas.body[Laplas.itek].arr_s_power,j+1);
         end;

         if (FormatSettings.DecimalSeparator='.') then
         begin
            // заменить все запятые в файле на точки.
            sub:=StringReplace(sub,',','.',[rfReplaceAll]);
         end;
         if (FormatSettings.DecimalSeparator=',') then
         begin
            // заменить все запятые в файле на точки.
            sub:=StringReplace(sub,'.',',',[rfReplaceAll]);
         end;
         if (identifier=1) then
         begin
            // thermal conductivity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]:=StrToFloat(sub);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]:=StrToFloat(sub);
         end;
         if (identifier=2) then
         begin
            // power
            Laplas.body[Laplas.itek].temp_power[j]:=StrToFloat(sub);
         end;

         s:=Copy(s,Pos(' ',s)+1,length(s));
         s:=Trim(s);
         if (FormatSettings.DecimalSeparator='.') then
         begin
            // заменить все запятые в файле на точки.
            s:=StringReplace(s,',','.',[rfReplaceAll]);
         end;
         if (FormatSettings.DecimalSeparator=',') then
         begin
            // заменить все запятые в файле на точки.
            s:=StringReplace(s,'.',',',[rfReplaceAll]);
         end;
         if (identifier=1) then
         begin
            // thermal conductivity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[j]:=StrToFloat(s);
         end;
         if (identifier=0) then
         begin
            // heat capacity
            Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[j]:=StrToFloat(s);
         end;
         if (identifier=2) then
         begin
            // power
            Laplas.body[Laplas.itek].arr_s_power[j]:=Trim(s);
            //Laplas.body[Laplas.itek].arr_power[j]:=StrToFloat(s);
            if (bOk) then
            begin
               Laplas.body[Laplas.itek].arr_power[j]:=FormVariables.my_real_convert(s,bOk);
            end;
         end;
         inc(j);
      end;
   end;


   if (identifier=1) then
   begin
      // thermal conductivity
      Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam:=j;
   end;
   if (identifier=0) then
   begin
      // heat capacity
      Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp:=j;
   end;
   if (identifier=2) then
   begin
      // power
      Laplas.body[Laplas.itek].n_power:=j;
   end;
    if (ComboBoxtemperatureUnit.ItemIndex=1) then
   begin
     // K 2 C
     if (identifier=1) then
     begin
        // thermal conductivity
        for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam-1 do
        begin
           Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j]-273.0;
        end;
     end;
     if (identifier=0) then
     begin
        // heat capacity
        for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp-1 do
        begin
           Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j]-273.0;
        end;
     end;
     if (identifier=2) then
     begin
        // power
        for j := 0 to Laplas.body[Laplas.itek].n_power-1 do
        begin
           Laplas.body[Laplas.itek].temp_power[j]:=Laplas.body[Laplas.itek].temp_power[j]-273.0;
        end;
     end;
   end;
   // Визуализация пользовательской зависимости.
   // TODO:
    frmRectangularPlot.cht1.Axes.Left.Logarithmic:=false;

    frmRectangularPlot.cht1.Title.Text[0]:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
    if (identifier=1) then
    begin
      // thermal conductivity
      frmRectangularPlot.cht1.LeftAxis.Title.Caption:='thermal conductivity, W/(m*K)';
    end;
     if (identifier=0) then
    begin
      // heat capacity
      frmRectangularPlot.cht1.LeftAxis.Title.Caption:='heat capacity, J/(kg*K)';
    end;
     if (identifier=2) then
    begin
      // temperature depend power
      frmRectangularPlot.cht1.LeftAxis.Title.Caption:='temperature depend power, W';
    end;
    frmRectangularPlot.cht1.BottomAxis.Title.Caption:='Temperature, C';
    frmRectangularPlot.cht1.Series[0].Clear;
    if (identifier=1) then
    begin
      // thermal conductivity
      for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam-1 do
      begin
         frmRectangularPlot.cht1.Series[0].AddXY(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[j],Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[j],'',clRed);
      end;
       if (Memopiecewiseproperties.Lines.Count=1) then
       begin
          frmRectangularPlot.cht1.Series[0].AddXY(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_lam[0]+100.0,Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[0],'',clRed);
       end;
    end;
    if (identifier=0) then
    begin
      // heat capacity
      for j := 0 to Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp-1 do
      begin
         frmRectangularPlot.cht1.Series[0].AddXY(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[j],Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[j],'',clRed);
      end;
      if (Memopiecewiseproperties.Lines.Count=1) then
       begin
           frmRectangularPlot.cht1.Series[0].AddXY(Laplas.workmat[Laplas.body[Laplas.itek].imatid].temp_cp[0]+100.0,Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[0],'',clRed);
       end;
    end;
    if (identifier=2) then
    begin
      // power as function temperature.
      for j := 0 to Laplas.body[Laplas.itek].n_power-1 do
      begin
         frmRectangularPlot.cht1.Series[0].AddXY(Laplas.body[Laplas.itek].temp_power[j],Laplas.body[Laplas.itek].arr_power[j],'',clRed);
      end;
       if (Memopiecewiseproperties.Lines.Count=1) then
       begin
          frmRectangularPlot.cht1.Series[0].AddXY(Laplas.body[Laplas.itek].temp_power[0]+100.0,Laplas.body[Laplas.itek].arr_power[0],'',clRed);
       end;
    end;
    frmRectangularPlot.ShowModal;
end;

end.
