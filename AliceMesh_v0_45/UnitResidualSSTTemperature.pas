unit UnitResidualSSTTemperature;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, VCLTee.TeEngine,
  VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidualSSTTemp = class(TForm)
    Chart1: TChart;
    Timer1: TTimer;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Series6: TFastLineSeries;
    Series7: TFastLineSeries;
    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure Timer1Timer(Sender: TObject);
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    brun_visibleSSTTemp : Boolean;
  end;

var
  FormResidualSSTTemp: TFormResidualSSTTemp;

implementation

{$R *.dfm}

uses VisualUnit, UnitEQGD;

// ������ ����� �������������.
procedure TFormResidualSSTTemp.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;

// ��������� �������� �����.
procedure TFormResidualSSTTemp.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidualSSTTemp.ClientHeight;
   Chart1.Width:=FormResidualSSTTemp.ClientWidth;
end;

procedure TFormResidualSSTTemp.Timer1Timer(Sender: TObject);
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
       if brun_visibleSSTTemp then
       begin

       if (Laplas.egddata.itemper>0) then
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
             FormResidualSSTTemp.Chart1.SeriesList[0].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[1].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[2].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[3].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[4].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[5].Clear;
             FormResidualSSTTemp.Chart1.SeriesList[6].Clear;

             if (f.Count>9) then
             begin
                istart:=6;
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
                FormResidualSSTTemp.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
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
                   FormResidualSSTTemp.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);

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
                      FormResidualSSTTemp.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
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
                         FormResidualSSTTemp.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);
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
                            FormResidualSSTTemp.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);

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
                               FormResidualSSTTemp.Chart1.SeriesList[5].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                                  FormResidualSSTTemp.Chart1.SeriesList[6].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                               end;
                            end
                             else
                            begin
                               // TODO
                               // ����� ������ ����� ������ ��� ��������.
                            end;
                         end;
                      end;
                  end;
               end
                else
               begin
                  // TODO
                  // ����� ������.
               end;
               //FormResidualSSTTemp.Chart1.LeftAxis.Minimum:=1e-3*fmin;
               //FormResidualSSTTemp.Chart1.LeftAxis.Maximum:=1e3*fmax;
               FormResidualSSTTemp.Chart1.LeftAxis.Minimum:=0.5*fmin;
               FormResidualSSTTemp.Chart1.LeftAxis.Maximum:=fmax;
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
        brun_visibleSSTTemp:=false;

      end;

    f.Free;
   end;
   end;
end;

end.
