unit Unitresidual2;
// ����������� ������� � ������ ����� ��� �������� � ��������� �������������.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, VCLTee.TeEngine,
  VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart, Vcl.AppEvnts;

type
  TFormresidual2 = class(TForm)
    cht1: TChart;
    Series3: TLineSeries;
    Series4: TLineSeries;
    Series5: TLineSeries;
    Timer1: TTimer;
    Series2: TLineSeries;
    Series1: TLineSeries;
    ApplicationEvents1: TApplicationEvents;
    procedure Timer1Timer(Sender: TObject);
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }

  public
    { Public declarations }
    brun_visible2 : Boolean;
  end;

var
  Formresidual2: TFormresidual2;

implementation

{$R *.dfm}

uses VisualUnit, UnitEQGD;

// ������ ����� �������������.
procedure TFormresidual2.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
     if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
 then
  msg.message:=0;
end;

// ��������� �������� �����.
procedure TFormresidual2.FormResize(Sender: TObject);
begin
   cht1.Height:=Formresidual2.ClientHeight;
   cht1.Width:=Formresidual2.ClientWidth;
end;

procedure TFormresidual2.Timer1Timer(Sender: TObject);
var
   f : TStringList; // ���������� ���� ������ TStringList
   i : Integer;
   fmin, fmax : Real;
   s, sub, subx : string;
   istart : Integer;

begin
     if (Laplas.ecology_btn) then
   begin
   if (EGDForm.ComboBoxTemperature.ItemIndex=1) then
   begin

    istart:=2;

    // �������� ����� ����������� ������ �������.
    f:=TStringList.Create();

      try

      if brun_visible2 then
      begin

        if (Laplas.egddata.itemper=1) then
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
         Formresidual2.cht1.SeriesList[0].Clear;
         Formresidual2.cht1.SeriesList[1].Clear;
         Formresidual2.cht1.SeriesList[2].Clear;
         Formresidual2.cht1.SeriesList[3].Clear;
         Formresidual2.cht1.SeriesList[4].Clear;

         if (f.Count>9) then
         begin
            istart:=7;
         end;

         for i:=istart to f.Count-1 do
         begin
            fmin:=20.0;
            fmax:=120.0;
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
            Formresidual2.cht1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
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
               Formresidual2.cht1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
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
                  Formresidual2.cht1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);
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
                     Formresidual2.cht1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                     s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                     sub:=s;
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
                         Formresidual2.cht1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clBlack);
                     end
                      else
                     begin
                        // TODO
                        // ����� ������ ����� ������ ��� ��������.
                     end;
                  end
                  else
                  begin
                      // TODO
                      // ����� ������ ����� ������ ��� ��������.
                  end;

                end
                else
                begin
                   // TODO
                   // ����� ������ ����� ���� ������ ��������.
                end;
            end
            else
            begin
               // TODO
               // ����� ������.
            end;
            //Formresidual2.cht1.LeftAxis.Minimum:=1e-3*fmin;
            //Formresidual2.cht1.LeftAxis.Maximum:=1e3*fmax;
            Formresidual2.cht1.LeftAxis.Minimum:=0.5*fmin;
            Formresidual2.cht1.LeftAxis.Maximum:=fmax;
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
        brun_visible2:=false;

      end;

    f.Free;
   end;
   end;
end;

end.
