unit UnitResidual_Langtry_Menter;
// Отображает невязки в ходе расчёта для модели Лангтрии и Ментора.
// При решении уравнений без температуры.
// 16.18   21.01.2021

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, VclTee.TeeGDIPlus, Vcl.ExtCtrls,
  VCLTee.TeEngine, VCLTee.Series, VCLTee.TeeProcs, VCLTee.Chart;

type
  TFormResidual_Langtry_Menter = class(TForm)
    Chart1: TChart;
    Series1: TFastLineSeries;
    Series2: TFastLineSeries;
    Series3: TFastLineSeries;
    Series4: TFastLineSeries;
    Series5: TFastLineSeries;
    Series6: TFastLineSeries;
    Series7: TFastLineSeries;
    Series8: TFastLineSeries;
    Timer1: TTimer;

    procedure ApplicationEvents1Message(var Msg: tagMSG; var Handled: Boolean);
    procedure Timer1Timer(Sender: TObject);
    procedure FormResize(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    brun_visible_Langtry_Menter : Boolean;
  end;

var
  FormResidual_Langtry_Menter: TFormResidual_Langtry_Menter;

implementation

{$R *.dfm}

uses VisualUnit, Math, UnitEQGD;

// Запрет форме сворачиваться.
procedure TFormResidual_Langtry_Menter.ApplicationEvents1Message(var Msg: tagMSG;
  var Handled: Boolean);
begin
    if (msg.message=WM_SYSCOMMAND)and(msg.wParam=SC_MINIMIZE)
   then
      msg.message:=0;
end;

// Изменение размеров формы.
procedure TFormResidual_Langtry_Menter.FormCreate(Sender: TObject);
begin
    // первые две строки нужно пропустить.
             FormResidual_Langtry_Menter.Chart1.SeriesList[0].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[1].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[2].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[3].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[4].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[5].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[6].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[7].Clear;
end;

procedure TFormResidual_Langtry_Menter.FormResize(Sender: TObject);
begin
   Chart1.Height:=FormResidual_Langtry_Menter.ClientHeight;
   Chart1.Width:=FormResidual_Langtry_Menter.ClientWidth;
end;

procedure TFormResidual_Langtry_Menter.Timer1Timer(Sender: TObject);
var
   f : TStringList; // переменная типа объект TStringList
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

    // Действие будет происходить каждую секунду.
    f:=TStringList.Create();

    try
       if brun_visible_Langtry_Menter then
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
             FormResidual_Langtry_Menter.Chart1.SeriesList[0].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[1].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[2].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[3].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[4].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[5].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[6].Clear;
             FormResidual_Langtry_Menter.Chart1.SeriesList[7].Clear;


            // FormResidual_Langtry_Menter.Chart1.LeftAxis.Minimum:=1.0e-14;
            // FormResidual_Langtry_Menter.Chart1.LeftAxis.Maximum:=1.0e5;

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
                FormResidual_Langtry_Menter.Chart1.SeriesList[0].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clred);
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
                   FormResidual_Langtry_Menter.Chart1.SeriesList[1].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);

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
                      FormResidual_Langtry_Menter.Chart1.SeriesList[2].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clgreen);
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
                            FormResidual_Langtry_Menter.Chart1.SeriesList[3].AddXY(StrToFloat(subx),1.0,subx,clblue);
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
                            FormResidual_Langtry_Menter.Chart1.SeriesList[3].AddXY(StrToFloat(subx),StrToFloat(sub)/m1,subx,clblue);
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
                            FormResidual_Langtry_Menter.Chart1.SeriesList[4].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                            FormResidual_Langtry_Menter.Chart1.SeriesList[5].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                            FormResidual_Langtry_Menter.Chart1.SeriesList[6].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);


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
                               FormResidual_Langtry_Menter.Chart1.SeriesList[7].AddXY(StrToFloat(subx),StrToFloat(sub),subx,clOlive);
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
                   FormResidual_Langtry_Menter.Chart1.LeftAxis.Minimum:=fmin-1.0e3;
                   FormResidual_Langtry_Menter.Chart1.LeftAxis.Maximum:=fmax+1.0e3;
                end
                else
                begin
                   FormResidual_Langtry_Menter.Chart1.LeftAxis.Minimum:=1e-3*fmin;
                   FormResidual_Langtry_Menter.Chart1.LeftAxis.Maximum:=1e3*fmax;
                end;
                end
                else
                begin
                   if (fmax<= 1.005*fmin) then
                   begin
                      FormResidual_Langtry_Menter.Chart1.LeftAxis.Minimum:=fmin-1.0e3;
                      FormResidual_Langtry_Menter.Chart1.LeftAxis.Maximum:=fmax+1.0e3;
                   end
                    else
                   begin
                      FormResidual_Langtry_Menter.Chart1.LeftAxis.Minimum:=0.5*fmin;
                      FormResidual_Langtry_Menter.Chart1.LeftAxis.Maximum:=fmax;
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
       brun_visible_Langtry_Menter:=false;

    end;

    f.Free;
   end;
   end;
end;


end.
