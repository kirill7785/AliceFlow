unit UnitresidualPlotSpallartAllmares;
// ���������� ������� � ���� ��������������� ��������
// ��� cfd ������ ��� ����������� �� � ������������
// �������������� ��������� �� ������ ��������-���������.
// �������� 2019 ����.


interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, VCLTee.TeEngine,
  VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidualSpallart_Allmares = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Timer1: TTimer;
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure Timer1Timer(Sender: TObject);
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    brun_visibleSA : Boolean;
  end;

var
  FormResidualSpallart_Allmares: TFormResidualSpallart_Allmares;

implementation

{$R *.dfm}

uses VisualUnit, UnitResidualSATemp2, UnitEQGD;

// ������ ����� �������������.
procedure TFormResidualSpallart_Allmares.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

// ��������� �������� �����.
procedure TFormResidualSpallart_Allmares.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidualSpallart_Allmares.ClientHeight;
   Chart1.Width:=FormResidualSpallart_Allmares.ClientWidth;
end;

procedure TFormResidualSpallart_Allmares.Timer1Timer(Sender: TObject);
var
   f : TStringList; // ���������� ���� ������ TStringList
   i : Integer;
   fmin, fmax, m1 : Real;
   s, sub, subx : string;
   istart : Integer;
begin
     if (Laplas.ecology_btn) then
   begin
   if (EGDForm.ComboBoxTemperature.ItemIndex=0) then
   begin
    m1:=1.0;
    istart:=2;

    // �������� ����� ����������� ������ �������.
    f:=TStringList.Create();

      try
       if brun_visibleSA then
       begin

       if (Laplas.egddata.itemper=0) then
     begin
     if (FileExists('statistic_convergence.txt')) then
      begin
         f.LoadFromFile('statistic_convergence.txt');

          if (FormatSettings.DecimalSeparator=',') then
          begin
             // �������� ��� ����� � ����� �� �������.
             for i:=0 to f.Count-1 do
             begin
                s:=f.Strings[i];
                f.Strings[i]:=StringReplace(s,'.',',',[rfReplaceAll]);
             end;
          end;

         // ������ ��� ������ ����� ����������.
         FormResidualSpallart_Allmares.Chart1.SeriesList[0].Clear;
         FormResidualSpallart_Allmares.Chart1.SeriesList[1].Clear;
         FormResidualSpallart_Allmares.Chart1.SeriesList[2].Clear;
         FormResidualSpallart_Allmares.Chart1.SeriesList[3].Clear;
         FormResidualSpallart_Allmares.Chart1.SeriesList[4].Clear;

         if (f.Count>9) then
         begin
            istart:=7;
         end;

         for i:=istart to f.Count-1 do
         begin
            fmin:=20.0;

            fmax:=1.2;
            if (i<7) then
            begin
               // ��������� ������� �������� ��������.
               fmax:=120.0;
            end;

            s:=Trim(f.Strings[i]);
            subx:=Trim(Copy(s,1,Pos(' ',s)));
            s:=Trim(Copy(s,Pos(' ',s),Length(s)));
            sub:=Trim(Copy(s,1,Pos(' ',s)));
            if (StrToFloat(sub)<fmin) then
            begin
               fmin:=StrToFloat(sub);
            end;
            if (StrToFloat(sub)>fmax) then
            begin
               fmax:=StrToFloat(sub);
            end;
            FormResidualSpallart_Allmares.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
            s:=Trim(Copy(s,Pos(' ',s),Length(s)));
            sub:=Trim(Copy(s,1,Pos(' ',s)));
            if (length(sub)>0) then
            begin
               if (StrToFloat(sub)<fmin) then
               begin
                  fmin:=StrToFloat(sub);
               end;
               if (StrToFloat(sub)>fmax) then
               begin
                  fmax:=StrToFloat(sub);
               end;
               FormResidualSpallart_Allmares.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
               s:=Trim(Copy(s,Pos(' ',s),Length(s)));
               sub:=Trim(Copy(s,1,Pos(' ',s)));
               if (length(sub)>0) then
               begin
                  if (StrToFloat(sub)<fmin) then
                  begin
                     fmin:=StrToFloat(sub);
                  end;
                  if (StrToFloat(sub)>fmax) then
                  begin
                     fmax:=StrToFloat(sub);
                  end;
                  FormResidualSpallart_Allmares.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);
                  s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                  if (pos(' ',Trim(s))=0) then
                  begin
                    sub:=s;
                  end
                  else
                  begin
                     sub:=Trim(Copy(s,1,Pos(' ',s)));
                  end;
                  if (length(sub)>0) then
                  begin
                     if (i=istart) then
                     begin
                        // ����������������� continity.
                        m1:=StrToFloat(sub);
                        if (m1<1.0e-20) then
                        begin
                           // ������ �� ������� �� ����.
                           m1:=1.0;
                        end;
                        FormResidualSpallart_Allmares.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub)/m1,subx,clOlive);
                     end
                     else
                     begin
                        if (StrToFloat(sub)/m1<fmin) then
                        begin
                           fmin:=StrToFloat(sub)/m1;
                        end;
                        if (StrToFloat(sub)/m1>fmax) then
                        begin
                           fmax:=StrToFloat(sub)/m1;
                        end;
                        FormResidualSpallart_Allmares.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub)/m1,subx,clOlive);
                     end;
                     s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                     if (pos(' ',Trim(s))=0) then
                     begin
                       sub:=s;
                     end
                     else
                     begin
                        sub:=Trim(Copy(s,1,Pos(' ',s)));
                     end;
                     if (length(sub)>0) then
                     begin
                        if (StrToFloat(sub)<fmin) then
                        begin
                           fmin:=StrToFloat(sub);
                        end;
                        if (StrToFloat(sub)>fmax) then
                        begin
                           fmax:=StrToFloat(sub);
                        end;
                        FormResidualSpallart_Allmares.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                     end;

                  end
                  else
                  begin
                      // TODO
                      // ����� ������ ����� ������ ��� ��������.
                  end;

                end;
            end
            else
            begin
               // TODO
               // ����� ������.
            end;
            if (f.Count<=9) then
            begin
               FormResidualSpallart_Allmares.Chart1.LeftAxis.Minimum:=1e-3*fmin;
               FormResidualSpallart_Allmares.Chart1.LeftAxis.Maximum:=1e3*fmax;
            end
            else
            begin
               FormResidualSpallart_Allmares.Chart1.LeftAxis.Minimum:=0.5*fmin;
               FormResidualSpallart_Allmares.Chart1.LeftAxis.Maximum:=fmax;
            end;
         end;
         // ��� ������� ��������� �����, ��� ����� ��� ����������
         // �� ��� ����� ��������� ��������� ����������.
         //Formresidual.Show;
      end
      else
      begin
         // ���� ���� �� ������ �� �� ���� ������ ��������� ���������.
         //MainMemo.Lines.Add('statistic_convergence.txt not found.');
      end;
     end;
       end;

     except
        brun_visibleSA:=false;

      end;

    f.Free;
   end;
   end;
end;

end.
