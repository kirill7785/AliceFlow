unit UnitResidual_Langtry_Menter_Temp;
// Отображает невязки в ходе расчёта для модели Лангтрии и Ментора.
// При решении уравнений с учётом температуры.
// 17.16   23.01.2021

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, Vcl.ExtCtrls,
  VCLTee.TeEngine, VCLTee.Series, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidual_Lagtry_Menter_Temp = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;  // X-Vel
    Series2: TFastLineSeries;  // Y-Vel
    Series3: TFastLineSeries;  // Z-Vel
    Series4: TFastLineSeries;  // continity
    Series5: TFastLineSeries;  // Температура
    Series6: TFastLineSeries;  // k
    Series7: TFastLineSeries;  // omega
    Series8: TFastLineSeries;  // gamma
    Series9: TFastLineSeries;  // ReTheta
    Timer1: TTimer;

    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure Timer1Timer(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure FormResize(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    brun_visible_Langtry_Menter_Temp : Boolean;
  end;

var
  FormResidual_Lagtry_Menter_Temp: TFormResidual_Lagtry_Menter_Temp;

implementation

{$R *.dfm}


uses VisualUnit, Math, UnitEQGD;

// Запрет форме сворачиваться.
procedure TFormResidual_Lagtry_Menter_Temp.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;

// Изменение размеров формы.
procedure TFormResidual_Lagtry_Menter_Temp.FormCreate(Sender: TObject);
begin
    // первые две строки нужно пропустить.
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[0].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[1].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[2].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[3].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[4].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[5].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[6].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[7].Clear;
    FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[8].Clear;
end;

procedure TFormResidual_Lagtry_Menter_Temp.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidual_Lagtry_Menter_Temp.ClientHeight;
   Chart1.Width:=FormResidual_Lagtry_Menter_Temp.ClientWidth;
end;

procedure TFormResidual_Lagtry_Menter_Temp.Timer1Timer(Sender: TObject);
var
   f : TStringList; // переменная типа объект TStringList
   i : Integer;
   fmin, fmax, m1 : Real;
   s, sub, subx : string;
   istart : Integer;

begin
 if (Laplas.ecology_btn) then
   begin
   if (EGDForm.ComboBoxTemperature.ItemIndex=1) then
   begin

   m1:=1.0;
   istart:=2;

    // Действие будет происходить каждую секунду.
    f:=TStringList.Create();

    try
       if brun_visible_Langtry_Menter_Temp then
       begin

        if (Laplas.egddata.itemper=0) then
        begin
          if (FileExists('statistic_convergence.txt')) then
          begin
             f.LoadFromFile('statistic_convergence.txt');

             if (FormatSettings.DecimalSeparator=',') then
             begin
                // заменить все точки в файле на запятые.
                for i:=0 to f.Count-1 do
                begin
                   s:=f.Strings[i];
                   f.Strings[i]:=StringReplace(s,'.',',',[rfReplaceAll]);
                end;
             end;

             // первые две строки нужно пропустить.
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[0].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[1].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[2].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[3].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[4].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[5].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[6].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[7].Clear;
             FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[8].Clear;


            // FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Minimum:=1.0e-14;
            // FormResidual_Lagtry_Menter_Temp.LeftAxis.Maximum:=1.0e5;

             if (f.Count>9) then
             begin
                istart:=7;
             end;

             for i:=istart to f.Count-1 do
             begin
                fmin:=20.0;
                fmax:=100.2;
                if (i<7) then
                begin
                   // Начальные сильные волнения неваязки.
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
                FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
                s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                sub:=Trim(Copy(s,1,Pos(' ',s)));
                if (length(sub)>0) then
                begin
                   if (StrToFloat(sub)<=fmin) then
                   begin
                      fmin:=Max(1.0e-12,StrToFloat(sub));
                   end;
                   if (StrToFloat(sub)>=fmax) then
                   begin
                      fmax:=StrToFloat(sub);
                   end;
                   FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);

                   s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                   sub:=Trim(Copy(s,1,Pos(' ',s)));
                   if (length(sub)>0) then
                   begin
                      if (StrToFloat(sub)<=fmin) then
                      begin
                         fmin:=Max(1.0e-12,StrToFloat(sub));
                      end;
                      if (StrToFloat(sub)>=fmax) then
                      begin
                         fmax:=StrToFloat(sub);
                      end;
                      FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
                      s:=Trim(Copy(s,Pos(' ',s),Length(s)));
                      sub:=Trim(Copy(s,1,Pos(' ',s)));
                      if (length(sub)>0) then
                      begin
                         if (i=istart) then
                         begin
                            // Ремасштабирование continity.
                            m1:=StrToFloat(sub);
                            if (m1<1.0e-12) then
                            begin
                               // защита от деления на ноль.
                               m1:=1.0;
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[3].AddXY(StrToFloat(subx),1.0,subx,clblue);
                         end
                         else
                         begin
                            if (StrToFloat(sub)/m1<=fmin) then
                            begin
                               fmin:=Max(1.0e-12,StrToFloat(sub)/m1);
                            end;
                            if (StrToFloat(sub)/m1>=fmax) then
                            begin
                               fmax:=StrToFloat(sub)/m1;
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub)/m1,subx,clblue);
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
                            if (StrToFloat(sub)<=fmin) then
                            begin
                               fmin:=Max(1.0e-12,StrToFloat(sub));
                            end;
                            if (StrToFloat(sub)>=fmax) then
                            begin
                               fmax:=StrToFloat(sub);
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                            if (StrToFloat(sub)<=fmin) then
                            begin
                               fmin:=Max(1.0e-12,StrToFloat(sub));
                            end;
                            if (StrToFloat(sub)>=fmax) then
                            begin
                               fmax:=StrToFloat(sub);
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[5].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                            if (StrToFloat(sub)<=fmin) then
                            begin
                               fmin:=Max(1.0e-12,StrToFloat(sub));
                            end;
                            if (StrToFloat(sub)>=fmax) then
                            begin
                               fmax:=StrToFloat(sub);
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[6].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                            if (StrToFloat(sub)<=fmin) then
                            begin
                               fmin:=Max(1.0e-12,StrToFloat(sub));
                            end;
                            if (StrToFloat(sub)>=fmax) then
                            begin
                               fmax:=StrToFloat(sub);
                            end;
                            FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[7].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);


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
                               if (StrToFloat(sub)<=fmin) then
                               begin
                                  fmin:=Max(1.0e-12,StrToFloat(sub));
                               end;
                               if (StrToFloat(sub)>=fmax) then
                               begin
                                  fmax:=StrToFloat(sub);
                               end;
                               FormResidual_Lagtry_Menter_Temp.Chart1.SeriesList[8].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
                            end;
                         end;
                      end;

                   end
                    else
                   begin
                      // TODO
                       // обрыв данных после первых трёх значений.
                   end;
                end;
             end;
          end
          else
          begin
             // TODO
             // обрыв данных.
          end;
                (*if (f.Count<=9) then
                begin
                if (fmax<= 1.005*fmin) then
                begin
                   FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Minimum:=fmin-1.0e3;
                   FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Maximum:=fmax+1.0e3;
                end
                else
                begin
                   FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Minimum:=1e-3*fmin;
                   FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Maximum:=1e3*fmax;
                end;
                end
                else
                begin
                   if (fmax<= 1.005*fmin) then
                   begin
                      FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Minimum:=fmin-1.0e3;
                      FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Maximum:=fmax+1.0e3;
                   end
                    else
                   begin
                      FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Minimum:=0.5*fmin;
                      FormResidual_Lagtry_Menter_Temp.Chart1.LeftAxis.Maximum:=fmax;
                   end;
                end; *)
        end;
             // Нам ненужно запускать форму, нам нужно при запущенной
             // из вне формы постоянно обновлять информацию.
             //Formresidual.Show;
       end
           else
          begin
             // Если файл не найден на не надо ничего считывать постоянно.
             //MainMemo.Lines.Add('statistic_convergence.txt not found.');
          end;



    except
       brun_visible_Langtry_Menter_Temp:=false;

    end;

    f.Free;
   end;
   end;

end;

end.
