unit UnitResidualMenterSST;
// форма отображени€ нев€зки во арем€ вычислени€.

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, VCLTee.TeEngine,
  VCLTee.Series, Vcl.ExtCtrls, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidualSST = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Series6: TFastLineSeries;
    Timer1: TTimer;

    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure Timer1Timer(Sender: TObject);
    procedure FormResize(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
    brun_visibleSST : Boolean;
  end;

var
  FormResidualSST: TFormResidualSST;

implementation

{$R *.dfm}

uses VisualUnit, UnitEQGD;

// «апрет форме сворачиватьс€.
procedure TFormResidualSST.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;


// »зменение размеров формы.
procedure TFormResidualSST.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidualSST.ClientHeight;
   Chart1.Width:=FormResidualSST.ClientWidth;
end;

procedure TFormResidualSST.Timer1Timer(Sender: TObject);
var
   f : TStringList; // переменна€ типа объект TStringList
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

    // ƒействие будет происходить каждую секунду.
    f:=TStringList.Create();

    try
       if brun_visibleSST then
       begin

        if (Laplas.egddata.itemper=0) then
        begin
          if (FileExists('statistic_convergence.txt')) then
          begin
             f.LoadFromFile('statistic_convergence.txt');

             if (FormatSettings.DecimalSeparator=',') then
             begin
                // заменить все точки в файле на зап€тые.
                for i:=0 to f.Count-1 do
                begin
                   s:=f.Strings[i];
                   f.Strings[i]:=StringReplace(s,'.',',',[rfReplaceAll]);
                end;
             end;

             // первые две строки нужно пропустить.
             FormResidualSST.Chart1.SeriesList[0].Clear;
             FormResidualSST.Chart1.SeriesList[1].Clear;
             FormResidualSST.Chart1.SeriesList[2].Clear;
             FormResidualSST.Chart1.SeriesList[3].Clear;
             FormResidualSST.Chart1.SeriesList[4].Clear;
             FormResidualSST.Chart1.SeriesList[5].Clear;

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
                   // Ќачальные сильные волнени€ нева€зки.
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
                FormResidualSST.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
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
                   FormResidualSST.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);

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
                      FormResidualSST.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
                      s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                      sub:=Trim(Copy(s,1,Pos(' ',s)));
                      if (length(sub)>0) then
                      begin
                         if (i=istart) then
                         begin
                            // –емасштабирование continity.
                            m1:=StrToFloat(sub);
                            if (m1<1.0e-20) then
                            begin
                               // защита от делени€ на ноль.
                               m1:=1.0;
                            end;
                            FormResidualSST.Chart1.SeriesList[3].AddXY(StrToFloat(subx),1.0,subx,clblue);
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
                            FormResidualSST.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub)/m1,subx,clblue);
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
                            FormResidualSST.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                               FormResidualSST.Chart1.SeriesList[5].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                            end;
                         end
                          else
                         begin
                            // TODO
                            // обрыв данных после первых трЄх значений.
                         end;
                      end;
                   end;
                end
                 else
                begin
                   // TODO
                   // обрыв данных.
                end;
                if (f.Count<=9) then
                begin
                   FormResidualSST.Chart1.LeftAxis.Minimum:=1e-3*fmin;
                   FormResidualSST.Chart1.LeftAxis.Maximum:=1e3*fmax;
                end
                else
                begin
                   FormResidualSST.Chart1.LeftAxis.Minimum:=0.5*fmin;
                   FormResidualSST.Chart1.LeftAxis.Maximum:=fmax;
                end;
             end;
             // Ќам ненужно запускать форму, нам нужно при запущенной
             // из вне формы посто€нно обновл€ть информацию.
             //Formresidual.Show;
          end
           else
          begin
             // ≈сли файл не найден на не надо ничего считывать посто€нно.
             //MainMemo.Lines.Add('statistic_convergence.txt not found.');
          end;
        end;
       end;

    except
       brun_visibleSST:=false;

    end;

    f.Free;
   end;
   end;
end;

end.
