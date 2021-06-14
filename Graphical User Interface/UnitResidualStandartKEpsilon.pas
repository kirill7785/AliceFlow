unit UnitResidualStandartKEpsilon;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, Vcl.ExtCtrls,
  VCLTee.TeEngine, VCLTee.Series, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidualStandartKEpsilon = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Series6: TFastLineSeries;
    Timer1: TTimer;
    procedure Timer1Timer(Sender: TObject);
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }
     procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
  public
    { Public declarations }
    brun_visibleKEpsilon : Boolean;
  end;

var
  FormResidualStandartKEpsilon: TFormResidualStandartKEpsilon;

implementation

{$R *.dfm}

uses VisualUnit, UnitEQGD;

// «апрет форме сворачиватьс€.
procedure TFormResidualStandartKEpsilon.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;

// изменение размеров формы.
procedure TFormResidualStandartKEpsilon.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidualStandartKEpsilon.Height;
   Chart1.Width:=FormResidualStandartKEpsilon.Width;
end;

procedure TFormResidualStandartKEpsilon.Timer1Timer(Sender: TObject);
var
   f : TStringList; // переменна€ типа объект TStringList
   i : Integer;
   fmin, fmax : Real;
   s, sub, subx : string;
begin
      if (Laplas.ecology_btn) then
   begin
   if (EGDForm.ComboBoxTemperature.ItemIndex=0) then
   begin
      // ƒействие будет происходить каждую секунду.
      f:=TStringList.Create();

      try
       if brun_visibleKEpsilon then
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
             FormResidualStandartKEpsilon.Chart1.SeriesList[0].Clear;
             FormResidualStandartKEpsilon.Chart1.SeriesList[1].Clear;
             FormResidualStandartKEpsilon.Chart1.SeriesList[2].Clear;
             FormResidualStandartKEpsilon.Chart1.SeriesList[3].Clear;
             FormResidualStandartKEpsilon.Chart1.SeriesList[4].Clear;
             FormResidualStandartKEpsilon.Chart1.SeriesList[5].Clear;
             for i:=2 to f.Count-1 do
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
                FormResidualStandartKEpsilon.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
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
                   FormResidualStandartKEpsilon.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);

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
                      FormResidualStandartKEpsilon.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
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
                         FormResidualStandartKEpsilon.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clblue);
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
                            FormResidualStandartKEpsilon.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                               FormResidualStandartKEpsilon.Chart1.SeriesList[5].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
            FormResidualStandartKEpsilon.Chart1.LeftAxis.Minimum:=1e-3*fmin;
            FormResidualStandartKEpsilon.Chart1.LeftAxis.Maximum:=1e3*fmax;
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
        brun_visibleKEpsilon:=false;

      end;

    f.Free;
   end;
   end;

end;

end.
