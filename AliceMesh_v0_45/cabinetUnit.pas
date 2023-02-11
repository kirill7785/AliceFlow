unit cabinetUnit;
// редактирование размеров кабинета

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TCabinetForm = class(TForm)
    PanelMain: TPanel;
    GroupBoxMCS: TGroupBox;
    GroupBoxangle: TGroupBox;
    Lalfa: TLabel;
    Lbeta: TLabel;
    EAlf: TEdit;
    EBet: TEdit;
    GroupBoxOrigin: TGroupBox;
    LabelxO: TLabel;
    LabelyO: TLabel;
    LabelzO: TLabel;
    EditXo: TEdit;
    EditYo: TEdit;
    EditZo: TEdit;
    GBCabinetMaterial: TGroupBox;
    RGCabinetMaterial: TRadioGroup;
    BEditMaterial: TButton;
    GBOperatingTemperature: TGroupBox;
    EOpTemp: TEdit;
    LCentiGrade: TLabel;
    lblAlphadiap: TLabel;
    lblBetadiap: TLabel;
    LGam: TLabel;
    EGam: TEdit;
    Label1: TLabel;
    Panel1: TPanel;
    LabelFilmCoefficient: TLabel;
    EditFilmCoefficient: TEdit;
    BCentiGrade: TButton;
    Label2: TLabel;
    ComboBoxFilmCoeff: TComboBox;
    Label3: TLabel;
    Panel2: TPanel;
    Label5: TLabel;
    ComboBoxGeometryTypeCabinet: TComboBox;
    GroupBoxcabinetsize: TGroupBox;
    LxS: TLabel;
    LyS: TLabel;
    LzS: TLabel;
    LxE: TLabel;
    LyE: TLabel;
    LzE: TLabel;
    Labelunit1: TLabel;
    Labelunit2: TLabel;
    Labelunit3: TLabel;
    Labelunit4: TLabel;
    Labelunit5: TLabel;
    Labelunit6: TLabel;
    ExS: TEdit;
    EyS: TEdit;
    EzS: TEdit;
    ExE: TEdit;
    EyE: TEdit;
    EzE: TEdit;
    LabelXoUnit: TLabel;
    LabelYoUnit: TLabel;
    LabelZoUnit: TLabel;
    procedure FormCreate(Sender: TObject);
    procedure RGCabinetMaterialClick(Sender: TObject);
    procedure BEditMaterialClick(Sender: TObject);
    procedure BCentiGradeClick(Sender: TObject);
    procedure ComboBoxFilmCoeffChange(Sender: TObject);
    procedure ComboBoxGeometryTypeCabinetChange(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  CabinetForm: TCabinetForm;

implementation
uses
     VisualUnit, UnitFlyuidMatLib, UnitUserDefinedFluidMaterial,
     UnitVariables, Unitdefineindentcabinet;
{$R *.dfm}




procedure TCabinetForm.FormCreate(Sender: TObject);
begin
   // вызывается при создании формы
   with Laplas.body[0] do
   begin
      ExS.Text:=sxS;
      EyS.Text:=syS;
      EzS.Text:=szS;
      ExE.Text:=sxE;
      EyE.Text:=syE;
      EzE.Text:=szE;
   end;
end;

// вводит координаты центра подвижной системы координат
procedure TCabinetForm.ComboBoxFilmCoeffChange(Sender: TObject);
   var
    s1, sforval : String;
    k : Integer;
    code : Integer;
    c : Real;
begin
    case ComboBoxFilmCoeff.ItemIndex of
    0 : begin
           // Адиабатическая стенка.
           EditFilmCoefficient.Visible:=false;
           LabelFilmCoefficient.Visible:=false;
           Label3.Visible:=false;
           Laplas.adiabatic_vs_heat_transfer_coeff:=0;
           Laplas.filmcoefficient:=0.0;
        end;
    1 : begin
           // Условие Ньютона Рихмана по умолчанию.
           EditFilmCoefficient.Visible:=true;
           LabelFilmCoefficient.Visible:=true;
           Label3.Visible:=true;
           Laplas.adiabatic_vs_heat_transfer_coeff:=1;
           s1:=Trim(EditFilmCoefficient.Text);
           sforval:='';
           sforval:=StringReplace(s1,',','.',[rfReplaceAll]);
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
              if (s1[k]=',') then sforval[k]:='.'
               else sforval[k]:=s1[k];
           end;
           EditFilmCoefficient.Text:=s1;
           val(sforval,c,code);
           if (code=0) then
            begin
               Laplas.filmcoefficient:=StrToFloat(EditFilmCoefficient.Text);
            end
             else
            begin
               EditFilmCoefficient.Text:=FloatToStr(Laplas.filmcoefficient);
               ShowMessage('Error! Please, input Film coefficient correctly...');
            end;
        end;
        2 : begin
               // Stefan Bolcman
               EditFilmCoefficient.Visible:=false;
               LabelFilmCoefficient.Visible:=false;
               Label3.Visible:=false;
               Laplas.adiabatic_vs_heat_transfer_coeff:=2;
               Laplas.filmcoefficient:=0.0;
            end;
           3 : begin
               // mixture condition
               EditFilmCoefficient.Visible:=true;
               LabelFilmCoefficient.Visible:=true;
               Label3.Visible:=true;
               Laplas.adiabatic_vs_heat_transfer_coeff:=3;
               s1:=Trim(EditFilmCoefficient.Text);
               sforval:='';
               sforval:=StringReplace(s1,',','.',[rfReplaceAll]);
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
                  if (s1[k]=',') then sforval[k]:='.'
                    else sforval[k]:=s1[k];
              end;
              EditFilmCoefficient.Text:=s1;
              val(sforval,c,code);
              if (code=0) then
              begin
                 Laplas.filmcoefficient:=StrToFloat(EditFilmCoefficient.Text);
              end
                else
              begin
                 EditFilmCoefficient.Text:=FloatToStr(Laplas.filmcoefficient);
                 ShowMessage('Error! Please, input Film coefficient correctly...');
              end;
               
        end;
    end;
end;


// Смена типа геометри кабинета.
procedure TCabinetForm.ComboBoxGeometryTypeCabinetChange(Sender: TObject);
var
    cab_geom, cab_geom2 : Visible_Line;
   lw_dec : Integer;
   i, j, i_1 : Integer; // счётчик цикла For.
   bmodelcheck_cab_hollow : Boolean; // для организации проверки корректности модели.
   bcontinue_delete : Boolean;  // для удаления всех стенок связанных с кабинетом.
   r : Real;
   bOk : Boolean;


// удаляет стенку с идентификатором id.
procedure delete_wall(id : Integer);
var
   il,jl : Integer;  // l - local
   wall_copy : array of TPlane;
begin
   if ((id >-1) and (id<Laplas.lw)) then
   begin
      SetLength(wall_copy,Laplas.lw-1);
      jl:=0;
      for il:=0 to (Laplas.lw-1) do
      begin
         if (il<>id) then
         begin
            wall_copy[jl]:=Laplas.wall[il];
            inc(jl);
         end;
      end;
      dec(Laplas.lw);
      SetLength(Laplas.wall,Laplas.lw);
      for il:=0 to (Laplas.lw-1) do
      begin
         Laplas.wall[il]:=wall_copy[il];
      end;
      SetLength(wall_copy,0); // освобождение памяти
      wall_copy:=nil;
   end;
end;

begin
    case ComboBoxGeometryTypeCabinet.ItemIndex of
       0 : begin
              // Fluid. Prism
              //  Изменить размер блока, при этом предложить
              // задать пользователю дистанцию от размеров объектов модели.
              // это нужно для задач обтекания на основе уравнений Навье-Стокса.
              // Вернуть из неактивности стенки пивязанные к размерам кабинета.
              if (Laplas.body[0].itype=2) then
              begin
                 // Ранее тип был HOLLOW, поэтому мы изменяем
                 // размер кабинета.
                 if (Laplas.lb=1) then
                 begin
                    // Мы только кабинет и имеем.
                    Laplas.body[0].xS:=-0.5;  Laplas.body[0].xE:=0.5;  // отображение параметризованной пользователем
                    Laplas.body[0].yS:=-0.5;  Laplas.body[0].yE:=0.5;  // геометрии
                    Laplas.body[0].zS:=-0.5;  Laplas.body[0].zE:=0.5;  // в реальные числа
                    if (FormatSettings.DecimalSeparator='.') then
                    begin
                       Laplas.body[0].sxS:='-0.5';  Laplas.body[0].sxE:='0.5';  // параметризованная
                       Laplas.body[0].syS:='-0.5';  Laplas.body[0].syE:='0.5';  // пользователем
                       Laplas.body[0].szS:='-0.5';  Laplas.body[0].szE:='0.5';  // геометрия
                    end;
                    if (FormatSettings.DecimalSeparator=',') then
                    begin
                       Laplas.body[0].sxS:='-0,5';  Laplas.body[0].sxE:='0,5';  // параметризованная
                       Laplas.body[0].syS:='-0,5';  Laplas.body[0].syE:='0,5';  // пользователем
                       Laplas.body[0].szS:='-0,5';  Laplas.body[0].szE:='0,5';  // геометрия
                    end;
                 end
                 else
                 begin
                    // Задание размеров кабинета с учтом остальных элементов
                    // присутствующих в модели.
                    // 1.Определить размеры элементов модели и если надо выдать
                    // предупреждение.
                    // 2.Вызвать окошко в котором указывается на какое расстояние надо
                    // отступить от размеров текущей модели, в том числе должно быть
                    // допустимо и нулевое расстояние.
                    // Показать пользователю наглядную картинку, чтобы он легче ориентировался.
                    // 3. Сразу предложить паттерн граничных условий на внешних стенках кабинета:
                    // 3.1. конвекция. 3.2. вынужденный обдув.
                    // HOLLOW BLOCK
                    bmodelcheck_cab_hollow:=true;

                    cab_geom.xS:=1.0e20;
                    cab_geom.xE:=-1.0e20;
                    cab_geom.yS:=1.0e20;
                    cab_geom.yE:=-1.0e20;
                    cab_geom.zS:=1.0e20;
                    cab_geom.zE:=-1.0e20;

                    for i:=1 to (Laplas.lb-1) do
                    begin
                       if (Laplas.body[i].bactivity) then
                       begin

                       with (Laplas.body[i]) do
                       begin
                          // 0 - PRISM; 1 - CYLINDER; 2 - POLYGON.
                          if (igeometry_type=2) then
                          begin
                             // 1 - XY, 2 - XZ, 3 - YZ.
                             case iPlane_obj2 of
                                1 : begin  //   1 - XY
                                    for i_1 := 0 to nsizei-1 do
                                    begin
                                       if (xi[i_1]<cab_geom.xS) then
                                       begin
                                          cab_geom.xS:=xi[i_1];
                                       end;
                                       if (yi[i_1]<cab_geom.yS) then
                                       begin
                                          cab_geom.yS:=yi[i_1];
                                       end;
                                       if (xi[i_1]>cab_geom.xE) then
                                       begin
                                          cab_geom.xE:=xi[i_1];
                                       end;
                                       if (yi[i_1]>cab_geom.yE) then
                                       begin
                                          cab_geom.yE:=yi[i_1];
                                       end;
                                       if (zi[i_1]<cab_geom.zS) then
                                       begin
                                          cab_geom.zS:=zi[i_1];
                                       end;
                                       if (zi[i_1]+hi[i_1]>cab_geom.zE) then
                                       begin
                                          cab_geom.zE:=zi[i_1]+hi[i_1];
                                       end;
                                    end;
                                end;
                                2 : begin // 2 - XZ
                                       if (xi[i_1]<cab_geom.xS) then
                                       begin
                                          cab_geom.xS:=xi[i_1];
                                       end;
                                       if (zi[i_1]<cab_geom.zS) then
                                       begin
                                          cab_geom.zS:=zi[i_1];
                                       end;
                                       if (xi[i_1]>cab_geom.xE) then
                                       begin
                                          cab_geom.xE:=xi[i_1];
                                       end;
                                       if (zi[i_1]>cab_geom.zE) then
                                       begin
                                          cab_geom.zE:=zi[i_1];
                                       end;
                                       if (yi[i_1]<cab_geom.yS) then
                                       begin
                                          cab_geom.yS:=yi[i_1];
                                       end;
                                       if (yi[i_1]+hi[i_1]>cab_geom.yE) then
                                       begin
                                          cab_geom.yE:=yi[i_1]+hi[i_1];
                                       end;
                                     end;
                                3 : begin  //  3 - YZ
                                        for i_1 := 0 to nsizei-1 do
                                        begin
                                           if (zi[i_1]<cab_geom.zS) then
                                           begin
                                              cab_geom.zS:=zi[i_1];
                                           end;
                                           if (yi[i_1]<cab_geom.yS) then
                                           begin
                                              cab_geom.yS:=yi[i_1];
                                           end;
                                           if (zi[i_1]>cab_geom.zE) then
                                           begin
                                              cab_geom.zE:=zi[i_1];
                                           end;
                                           if (yi[i_1]>cab_geom.yE) then
                                           begin
                                              cab_geom.yE:=yi[i_1];
                                           end;
                                           if (xi[i_1]<cab_geom.xS) then
                                           begin
                                              cab_geom.xS:=xi[i_1];
                                           end;
                                           if (xi[i_1]+hi[i_1]>cab_geom.xE) then
                                           begin
                                              cab_geom.xE:=xi[i_1]+hi[i_1];
                                           end;
                                        end;
                                    end;
                             end;
                          end
                          else
                          begin
                             if (xS<cab_geom.xS) then
                             begin
                                cab_geom.xS:=xS;
                             end;
                             if (yS<cab_geom.yS) then
                             begin
                                cab_geom.yS:=yS;
                             end;
                             if (zS<cab_geom.zS) then
                             begin
                                cab_geom.zS:=zS;
                             end;
                             if (xE>cab_geom.xE) then
                             begin
                                cab_geom.xE:=xE;
                             end;
                             if (yE>cab_geom.yE) then
                             begin
                                cab_geom.yE:=yE;
                             end;
                             if (zE>cab_geom.zE) then
                             begin
                                cab_geom.zE:=zE;
                             end;
                          end;
                       end;
                     end;
                    end;

                    cab_geom2.xS:=1.0e20;
                    cab_geom2.xE:=-1.0e20;
                    cab_geom2.yS:=1.0e20;
                    cab_geom2.yE:=-1.0e20;
                    cab_geom2.zS:=1.0e20;
                    cab_geom2.zE:=-1.0e20;
                    for i:=0 to (Laplas.ls-1) do
                    begin
                       with (Laplas.source[i]) do
                       begin
                          if (xS<cab_geom2.xS) then
                          begin
                             cab_geom2.xS:=xS;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' xS position is incorrect.');
                          end;
                          if (yS<cab_geom2.yS) then
                          begin
                             cab_geom2.yS:=yS;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' yS position is incorrect.');
                          end;
                          if (zS<cab_geom2.zS) then
                          begin
                             cab_geom2.zS:=zS;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' zS position is incorrect.');
                          end;
                          if (xE>cab_geom2.xE) then
                          begin
                             cab_geom2.xE:=xE;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' xE position is incorrect.');
                          end;
                          if (yE>cab_geom2.yE) then
                          begin
                             cab_geom2.yE:=yE;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' yE position is incorrect.');
                          end;
                          if (zE>cab_geom2.zE) then
                          begin
                             cab_geom2.zE:=zE;
                             Laplas.MainMemo.Lines.Add('error:'+Laplas.source[i].name+' zE position is incorrect.');
                          end;
                       end;
                     end;

                     for i:=0 to (Laplas.lw-1) do
                     begin
                        with (Laplas.wall[i]) do
                        begin
                           if (cabinet_depend=0) then
                           begin
                              if (xS<cab_geom2.xS) then
                              begin
                                 cab_geom2.xS:=xS;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' xS position is incorrect.');
                              end;
                              if (yS<cab_geom2.yS) then
                              begin
                                 cab_geom2.yS:=yS;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' yS position is incorrect.');
                              end;
                              if (zS<cab_geom2.zS) then
                              begin
                                 cab_geom2.zS:=zS;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' zS position is incorrect.');
                              end;
                              if (xE>cab_geom2.xE) then
                              begin
                                 cab_geom2.xE:=xE;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' xE position is incorrect.');
                              end;
                              if (yE>cab_geom2.yE) then
                              begin
                                 cab_geom2.yE:=yE;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' yE position is incorrect.');
                              end;
                              if (zE>cab_geom2.zE) then
                              begin
                                 cab_geom2.zE:=zE;
                                 Laplas.MainMemo.Lines.Add('error:'+Laplas.wall[i].name+' zE position is incorrect.');
                              end;
                           end;
                        end;
                     end;

                     if (cab_geom2.xS<cab_geom.xS) then
                     begin
                        bmodelcheck_cab_hollow:=false;
                     end;
                     if (cab_geom2.yS<cab_geom.yS) then
                     begin
                         bmodelcheck_cab_hollow:=false;
                     end;
                     if (cab_geom2.zS<cab_geom.zS) then
                     begin
                        bmodelcheck_cab_hollow:=false;
                     end;
                     if (cab_geom2.xE>cab_geom.xE) then
                     begin
                        bmodelcheck_cab_hollow:=false;
                     end;
                     if (cab_geom2.yE>cab_geom.yE) then
                     begin
                        bmodelcheck_cab_hollow:=false;
                     end;
                     if (cab_geom2.zE>cab_geom.zE) then
                     begin
                        bmodelcheck_cab_hollow:=false;
                     end;
                     if (not(bmodelcheck_cab_hollow)) then
                     begin
                        ShowMessage('Error. Your geometry model is incorrect. Dont Start Solver.');
                        Laplas.MainMemo.Lines.Add('Error. Your geometry model is incorrect. Dont Start Solver.');
                     end;
                     // Здесь нужно вызывать форму и спрашивать величину отступа.
                     Formcabinetindent.RadioGroup1.ItemIndex:=0; // без отступа по умолчанию.
                     Formcabinetindent.Edithonly.Text:=FloatToStr(0.0);

                     Formcabinetindent.Edithx.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edithy.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edithz.Text:=FloatToStr(0.0);

                     Formcabinetindent.Edit1.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edit2.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edit3.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edit4.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edit5.Text:=FloatToStr(0.0);
                     Formcabinetindent.Edit6.Text:=FloatToStr(0.0);

                     if (Laplas.ComboBoxlength.ItemIndex=1) then
                     begin
                        Formcabinetindent.Labelunion.Caption:='mm';

                        Formcabinetindent.Labelunion1.Caption:='mm';
                        Formcabinetindent.Labelunion2.Caption:='mm';
                        Formcabinetindent.Labelunion3.Caption:='mm';

                        Formcabinetindent.Label11.Caption:='mm';
                        Formcabinetindent.Label12.Caption:='mm';
                        Formcabinetindent.Label13.Caption:='mm';
                        Formcabinetindent.Label14.Caption:='mm';
                        Formcabinetindent.Label15.Caption:='mm';
                        Formcabinetindent.Label16.Caption:='mm';

                     end;
                     if (Laplas.ComboBoxlength.ItemIndex=0) then
                     begin
                        Formcabinetindent.Labelunion.Caption:='m';

                        Formcabinetindent.Labelunion1.Caption:='m';
                        Formcabinetindent.Labelunion2.Caption:='m';
                        Formcabinetindent.Labelunion3.Caption:='m';

                        Formcabinetindent.Label11.Caption:='m';
                        Formcabinetindent.Label12.Caption:='m';
                        Formcabinetindent.Label13.Caption:='m';
                        Formcabinetindent.Label14.Caption:='m';
                        Formcabinetindent.Label15.Caption:='m';
                        Formcabinetindent.Label16.Caption:='m';

                     end;
                     if (Laplas.ComboBoxlength.ItemIndex=2) then
                     begin
                        Formcabinetindent.Labelunion.Caption:='um';

                        Formcabinetindent.Labelunion1.Caption:='um';
                        Formcabinetindent.Labelunion2.Caption:='um';
                        Formcabinetindent.Labelunion3.Caption:='um';

                        Formcabinetindent.Label11.Caption:='um';
                        Formcabinetindent.Label12.Caption:='um';
                        Formcabinetindent.Label13.Caption:='um';
                        Formcabinetindent.Label14.Caption:='um';
                        Formcabinetindent.Label15.Caption:='um';
                        Formcabinetindent.Label16.Caption:='um';

                     end;
                     Formcabinetindent.ShowModal;
                     // Для случая нулевого отступа.
                     if (Formcabinetindent.RadioGroup1.ItemIndex=0) then
                     begin
                        Laplas.body[0].xS:=cab_geom.xS;  Laplas.body[0].xE:=cab_geom.xE;  // отображение параметризованной пользователем
                        Laplas.body[0].yS:=cab_geom.yS;  Laplas.body[0].yE:=cab_geom.yE;  // геометрии
                        Laplas.body[0].zS:=cab_geom.zS;  Laplas.body[0].zE:=cab_geom.zE;  // в реальные числа
                     end;
                     if (Formcabinetindent.RadioGroup1.ItemIndex=1) then
                     begin
                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edithonly.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edithonly.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;
                        Laplas.body[0].xS:=cab_geom.xS-r;  Laplas.body[0].xE:=cab_geom.xE+r;  // отображение параметризованной пользователем
                        Laplas.body[0].yS:=cab_geom.yS-r;  Laplas.body[0].yE:=cab_geom.yE+r;  // геометрии
                        Laplas.body[0].zS:=cab_geom.zS-r;  Laplas.body[0].zE:=cab_geom.zE+r;  // в реальные числа
                     end;
                     if (Formcabinetindent.RadioGroup1.ItemIndex=2) then
                     begin
                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edithx.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edithx.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;
                        Laplas.body[0].xS:=cab_geom.xS-r;  Laplas.body[0].xE:=cab_geom.xE+r;  // отображение параметризованной пользователем

                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edithy.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edithy.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                        Laplas.body[0].yS:=cab_geom.yS-r;  Laplas.body[0].yE:=cab_geom.yE+r;  // геометрии

                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edithz.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edithz.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                        Laplas.body[0].zS:=cab_geom.zS-r;  Laplas.body[0].zE:=cab_geom.zE+r;  // в реальные числа
                     end;

                     if (Formcabinetindent.RadioGroup1.ItemIndex=3) then
                     begin
                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit1.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit1.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;
                        Laplas.body[0].xS:=cab_geom.xS-r;

                         bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit2.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit2.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                         Laplas.body[0].xE:=cab_geom.xE+r;  // отображение параметризованной пользователем

                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit3.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit3.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                        Laplas.body[0].yS:=cab_geom.yS-r;

                          bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit4.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit4.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;


                         Laplas.body[0].yE:=cab_geom.yE+r;  // геометрии

                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit5.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit5.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                        Laplas.body[0].zS:=cab_geom.zS-r;

                        bOk:=true;
                        //r:=StrToFloat(Formcabinetindent.Edit6.Text);
                        r:=FormVariables.my_real_convert(Trim(Formcabinetindent.Edit6.Text),bOk);
                        if (not(bOk)) then
                        begin
                           r:=0;
                        end;

                        Laplas.body[0].zE:=cab_geom.zE+r;  // в реальные числа
                     end;


                     // Текстовый вид координат.
                     Laplas.body[0].sxS:=FloatToStr(Laplas.body[0].xS);
                     Laplas.body[0].sxE:=FloatToStr(Laplas.body[0].xE);  // отображение параметризованной пользователем
                     Laplas.body[0].syS:=FloatToStr(Laplas.body[0].yS);
                     Laplas.body[0].syE:=FloatToStr(Laplas.body[0].yE);  // геометрии
                     Laplas.body[0].szS:=FloatToStr(Laplas.body[0].zS);
                     Laplas.body[0].szE:=FloatToStr(Laplas.body[0].zE);  // в реальные числа


                     // Здесь также сразу можно определиться с гидродинамической задачей.
                     // Например сразу можно всё задать для расчёта естественной конвекции или
                     // обтекания.
                     // TODO.

                 end;
              end;
               // вызов редактирования свойств кабинета
              ExS.Text:=Laplas.body[0].sxS;   // параметризованные
              EyS.Text:=Laplas.body[0].syS;   // пользователем
              EzS.Text:=Laplas.body[0].szS;   // геометрические
              ExE.Text:=Laplas.body[0].sxE;   // размеры.
              EyE.Text:=Laplas.body[0].syE;
              EzE.Text:=Laplas.body[0].szE;
              GroupBoxCabinetSize.Visible:=true;
              GBcabinetmaterial.Visible:=true;
              Laplas.body[0].itype:=3; // Fluid
       end;
       1 : begin
              // HOLLOW
              // удалить все стенки из дерева элементов и из списка объектов которые привязаны
              // к размерам кабинета.

              bcontinue_delete:=true;
              while (bcontinue_delete) do
              begin
                 bcontinue_delete:=false;
                  for i:=0 to (Laplas.lw-1) do
                     begin
                        with (Laplas.wall[i]) do
                        begin
                           if (cabinet_depend>0) then
                           begin
                               // удалить стенку с номером i.
                               // в дереве тоже надо почистить.
                                // Чистка дерева.
                                for j := 0 to Laplas.MainTreeView.Items.Count-1 do
                                begin
                                   if(( Pos(Laplas.MainTreeView.Items[j].Text,Laplas.wall[i].name)=1)and(length(Laplas.MainTreeView.Items[j].Text)=length(Laplas.wall[i].name))) then
                                   begin
                                       Laplas.MainTreeView.Items[j].Selected:=true;
                                   end;
                                end;
                                if ((not((length(Laplas.MainTreeView.Selected.Text)=Length(Laplas.body[0].name))and(Pos(Laplas.MainTreeView.Selected.Text,Laplas.body[0].name)=1)))) then
                                begin
                                   Laplas.MainTreeView.Items.Delete(Laplas.MainTreeView.Selected);
                                end;
                               bcontinue_delete:=true;
                               Laplas.MainMemo.Lines.Add('mesage: wall '+Laplas.wall[i].name+' successfully removed.');
                               delete_wall(i);
                               break;
                           end;
                        end;
                     end;
              end;
              Laplas.cab_bound_condition.bminX:=false;
              Laplas.cab_bound_condition.bminY:=false;
              Laplas.cab_bound_condition.bminZ:=false;
              Laplas.cab_bound_condition.bmaxX:=false;
              Laplas.cab_bound_condition.bmaxY:=false;
              Laplas.cab_bound_condition.bmaxZ:=false;


              GroupBoxCabinetSize.Visible:=false;
              GBcabinetmaterial.Visible:=false;
              Laplas.body[0].itype:=2; // Hollow
       end;
    end;
end;

// ввод значений углов в радианах
procedure TCabinetForm.RGCabinetMaterialClick(Sender: TObject);
begin
   Application.MessageBox('Change happens only when you click Edit','Attantion!')
end;

// Редактирование свойств материала.
procedure TCabinetForm.BEditMaterialClick(Sender: TObject);
begin
   Laplas.itek:=0; // Cabinet
   if (RGCabinetMaterial.ItemIndex=0) then
      begin
         // FLUID Program Library
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            // инициализация
            FormFluidLibMat.ComboBoxFluidLibMaterial.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident-1;
             case FormFluidLibMat.ComboBoxFluidLibMaterial.ItemIndex of
              0 : begin
                     // Dry Air
                     FormFluidLibMat.LabelDensityValue.Caption:='1.10174';
                  end;
              1 : begin
                     // Water Liquid
                     FormFluidLibMat.LabelDensityValue.Caption:='874.525';
                  end;
            end;
         end;
         FormFluidLibMat.ShowModal; // библиотечные материалы
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=1) then
         begin
            case Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilibident of
              1 : begin
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Dry_Air';
                  end;
              2 : begin
                    Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat:='Water_Liquid';
                  end;
              end;
         end;
      end;

      // User-Defined Fluid material
      if (RGCabinetMaterial.ItemIndex=1) then
      begin
         FormUserDefinedFluidMaterial.CBTipPattern.ItemIndex:=0; // no pattern
          // FLUID User - Defined
         if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].blibmat=0) then
         begin
            // инициализация:
            FormUserDefinedFluidMaterial.CBRho.ItemIndex:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].bBoussinesq;
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].bBoussinesq=1) then
            begin
                // приближение Обербека Буссинеска.
                FormUserDefinedFluidMaterial.GBVolexpans.Visible:=true;
            end
            else
            begin
               FormUserDefinedFluidMaterial.GBVolexpans.Visible:=false;
            end;

            FormUserDefinedFluidMaterial.EMatName.Text:=Laplas.workmat[Laplas.body[Laplas.itek].imatid].namemat;
            FormUserDefinedFluidMaterial.ERho.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].rho);
            //FormUserDefinedFluidMaterial.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].cp);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_cp=1) then
            begin
               // Constant properties
               FormUserDefinedFluidMaterial.CBCp.ItemIndex:=0;
               FormUserDefinedFluidMaterial.CBCpChange(Sender);
               FormUserDefinedFluidMaterial.ECp.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_cp[0]);
            end
             else
            begin
               // Piecewise linear properties.
               FormUserDefinedFluidMaterial.CBCp.ItemIndex:=1;
               FormUserDefinedFluidMaterial.CBCpChange(Sender);
            end;
            //FormUserDefinedFluidMaterial.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].lambda);
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].n_lam=1) then
            begin
               // Constant properties.
               FormUserDefinedFluidMaterial.CBLam.ItemIndex:=0;
               FormUserDefinedFluidMaterial.CBLamChange(Sender);
               FormUserDefinedFluidMaterial.ELam.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].arr_lam[0]);
            end
             else
            begin
               // piecewise properties
              FormUserDefinedFluidMaterial.CBLam.ItemIndex:=1;
              FormUserDefinedFluidMaterial.CBLamChange(Sender);
            end;
            // динамическая вязкость
            if (Laplas.workmat[Laplas.body[Laplas.itek].imatid].ilawmu=0) then
            begin
               FormUserDefinedFluidMaterial.CBMu.ItemIndex:=0;
               FormUserDefinedFluidMaterial.LSIMu.Visible:=true;
               FormUserDefinedFluidMaterial.EMu.Visible:=true;
               FormUserDefinedFluidMaterial.EMu.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].mu);
            end
            else
            begin
               FormUserDefinedFluidMaterial.CBMu.ItemIndex:=1;
               FormUserDefinedFluidMaterial.EMu.Visible:=false;
               FormUserDefinedFluidMaterial.LSIMu.Visible:=false;
               FormUserDefinedFluidMaterial.EMu.Text:=' ';
            end;

            FormUserDefinedFluidMaterial.EBeta_T.Text:=FloatToStr(Laplas.workmat[Laplas.body[Laplas.itek].imatid].beta_t);
         end;
         FormUserDefinedFluidMaterial.ShowModal;
      end;
end;

// Ввод Operating Temperature
// Опорное значение температуры одно для всех расчётных областей.
procedure TCabinetForm.BCentiGradeClick(Sender: TObject);
var
   bOk, bOk2 : Boolean;
   code : Integer;
   c : Real;
   sforval : String;
   k : Integer;
   iAlf, iBet, iGam : Integer; // Целочисленные углы Эйлера.
   s1, s2, s3, s4, s5, s6 : String;
   r1, r2, r3, r4, r5, r6 : Real;
   i, itekcab23 : Integer;

begin

   s1:=Trim(EOpTemp.Text);
   sforval:='';
   sforval:=StringReplace(s1,',','.',[rfReplaceAll]);
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
      if (s1[k]=',') then sforval[k]:='.'
      else sforval[k]:=s1[k];
   end;
   EOpTemp.Text:=s1;

   bOk:=False;
  // val(EOpTemp.Text,c,code); // пытаемся получить число
   val(sforval,c,code);
   if (code=0) then
   begin
      Laplas.operatingtemperature:=StrToFloat(EOpTemp.Text);
      bOk:=true;
   end
    else
   begin
      EOpTemp.Text:=FloatToStr(Laplas.operatingtemperature);
      ShowMessage('Error! Please, input Temperature Ambient correctly...');
      bOk:=false;
   end;

    if (ComboBoxFilmCoeff.ItemIndex=1) then
    begin
       Laplas.adiabatic_vs_heat_transfer_coeff:=1;
       s1:=Trim(EditFilmCoefficient.Text);
       sforval:='';
       sforval:=StringReplace(s1,',','.',[rfReplaceAll]);
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
          if (s1[k]=',') then sforval[k]:='.'
          else sforval[k]:=s1[k];
       end;
       EditFilmCoefficient.Text:=s1;
       val(sforval,c,code);
       if (code=0) then
       begin
          Laplas.filmcoefficient:=StrToFloat(EditFilmCoefficient.Text);
       end
        else
       begin
          EditFilmCoefficient.Text:=FloatToStr(Laplas.filmcoefficient);
          ShowMessage('Error! Please, input Film coefficient correctly...');
          bOk:=false;
       end;
    end
    else if (ComboBoxFilmCoeff.ItemIndex=2) then
    begin
       Laplas.adiabatic_vs_heat_transfer_coeff:=2;
    end
    else if (ComboBoxFilmCoeff.ItemIndex=3) then
    begin
       Laplas.adiabatic_vs_heat_transfer_coeff:=3;
        s1:=Trim(EditFilmCoefficient.Text);
       sforval:='';
       sforval:=StringReplace(s1,',','.',[rfReplaceAll]);
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
          if (s1[k]=',') then sforval[k]:='.'
          else sforval[k]:=s1[k];
       end;
       EditFilmCoefficient.Text:=s1;
       val(sforval,c,code);
       if (code=0) then
       begin
          Laplas.filmcoefficient:=StrToFloat(EditFilmCoefficient.Text);
       end
        else
       begin
          EditFilmCoefficient.Text:=FloatToStr(Laplas.filmcoefficient);
          ShowMessage('Error! Please, input Film coefficient correctly...');
          bOk:=false;
       end;

    end
    else
    begin
       // Адиабатическая стенка.
       Laplas.adiabatic_vs_heat_transfer_coeff:=0;
    end;


    // Alpha, Beta, Gamma
     s1:=Trim(EAlf.Text);
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
   EAlf.Text:=s1;

    s1:=Trim(EBet.Text);
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
   EBet.Text:=s1;

    s1:=Trim(EGam.Text);
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
   EGam.Text:=s1;


   // ввод значений углов (углы вводятся в радианах !
   Laplas.Alf:=StrToFloat(EAlf.Text);
   Laplas.Bet:=StrToFloat(EBet.Text);
   Laplas.Gam:=StrToFloat(EGam.Text);

   iAlf:=Round((Laplas.Alf/3.141)*180);
   // конструкция приведённая ниже должна обеспечивать непрерывное
   // плавное вращение во всём спектре.
   if (iAlf>360) then
   begin
     iAlf:=(iAlf-360);
   end
   else
   if (iAlf<-360) then
   begin
     iAlf:=(iAlf+360);
   end;
   Laplas.Alf:=(iAlf/180)*3.141;

   iBet:=Round((Laplas.Bet/3.141)*180);
   // конструкция приведённая ниже должна обеспечивать непрерывное
   // плавное вращение во всём спектре.
    if (iBet>360) then
   begin
     iBet:=(iBet-360);
   end
   else
   if (iBet<-360) then
   begin
     iBet:=(iBet+360);
   end;

   Laplas.Bet:=(iBet/180)*3.141;

   iGam:=Round((Laplas.Gam/3.141)*180);
   // конструкция приведённая ниже должна обеспечивать непрерывное
   // плавное вращение во всём спектре.
    if (iGam>360) then
   begin
     iGam:=(iGam-360);
   end
   else
   if (iGam<-360) then
   begin
     iGam:=(iGam+360);
   end;

   Laplas.Gam:=(iGam/180)*3.141;


   EAlf.Text:=FloatToStr(Laplas.Alf);
   EBet.Text:=FloatToStr(Laplas.Bet);
   EGam.Text:=FloatToStr(Laplas.Gam);

   // мгновенная перерисовка.

   // 3
   // Центр системы координат.
    s1:=Trim(EditXo.Text);
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
   EditXo.Text:=s1;

    s1:=Trim(EditYo.Text);
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
   EditYo.Text:=s1;

    s1:=Trim(EditZo.Text);
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
   EditZo.Text:=s1;


   bOk2:=false;
  // val(EditXo.Text,c,code); // пытаемся получить число
  sforval:='';
  sforval:=StringReplace(EditXo.Text,',','.',[rfReplaceAll]);
   val(sforval,c,code);
   if (code=0) then
   begin
      //val(EditYo.Text,c,code); // пытаемся получить число
      sforval:='';
      sforval:=StringReplace(EditYo.Text,',','.',[rfReplaceAll]);
      val(sforval,c,code);
      if (code=0) then
      begin
         sforval:='';
         sforval:=StringReplace(EditZo.Text,',','.',[rfReplaceAll]);
         val(sforval,c,code);
         //val(EditZo.Text,c,code); // пытаемся получить число
         if (code=0) then
         begin
            bOk2:=true;
         end;
      end;
   end;

   if (bOk2) then
   begin
      // позиция центра подвижной системы координат
      Laplas.Oxc:=StrToFloat(EditXo.Text);
      Laplas.Oyc:=StrToFloat(EditYo.Text);
      Laplas.Ozc:=StrToFloat(EditZo.Text);
      // прорисовка после присваивания координат центральной точки
   end
    else
   begin
      EditXo.Text:=FloatToStr(Laplas.Oxc);
      EditYo.Text:=FloatToStr(Laplas.Oyc);
      EditZo.Text:=FloatToStr(Laplas.Ozc);
      ShowMessage('Error! Please, input Origin correctly...');
      bOk:=false;
   end;

   // 4
   // Размеры кабинета.
   // измененние свойств кабинета
    // инициализация :
   r1:=0.0;
   r2:=0.0;
   r3:=0.0;
   r4:=0.0;
   r5:=0.0;
   r6:=0.0;

   bOk2:=true; // признак правильности ввода

   // координаты блока

   s1:=ExS.Text;  // параметризованные
   s2:=EyS.Text;  // геометрические
   s3:=ExE.Text;  // размеры
   s4:=EyE.Text;  // заданные
   s5:=EzS.Text;  // пользователем
   s6:=EzE.Text;

   if (FormatSettings.DecimalSeparator='.') then
   begin
      for i:=1 to length(s1) do
      begin
         if (s1[i]=',') then s1[i]:='.';
      end;
      Exs.Text:=Trim(s1);
      for i:=1 to length(s2) do
      begin
         if (s2[i]=',') then s2[i]:='.';
      end;
      Eys.Text:=Trim(s2);
      for i:=1 to length(s3) do
      begin
         if (s3[i]=',') then s3[i]:='.';
      end;
      ExE.Text:=Trim(s3);
      for i:=1 to length(s4) do
      begin
         if (s4[i]=',') then s4[i]:='.';
      end;
      EyE.Text:=Trim(s4);
      for i:=1 to length(s5) do
      begin
         if (s5[i]=',') then s5[i]:='.';
      end;
      EzS.Text:=Trim(s5);
      for i:=1 to length(s6) do
      begin
         if (s6[i]=',') then s6[i]:='.';
      end;
      EzE.Text:=Trim(s6);
   end;

    if (FormatSettings.DecimalSeparator=',') then
   begin
      for i:=1 to length(s1) do
      begin
         if (s1[i]='.') then s1[i]:=',';
      end;
      Exs.Text:=Trim(s1);
      for i:=1 to length(s2) do
      begin
         if (s2[i]='.') then s2[i]:=',';
      end;
      Eys.Text:=Trim(s2);
      for i:=1 to length(s3) do
      begin
         if (s3[i]='.') then s3[i]:=',';
      end;
      ExE.Text:=Trim(s3);
      for i:=1 to length(s4) do
      begin
         if (s4[i]='.') then s4[i]:=',';
      end;
      EyE.Text:=Trim(s4);
      for i:=1 to length(s5) do
      begin
         if (s5[i]='.') then s5[i]:=',';
      end;
      EzS.Text:=Trim(s5);
      for i:=1 to length(s6) do
      begin
         if (s6[i]='.') then s6[i]:=',';
      end;
      EzE.Text:=Trim(s6);
   end;


   if bOk2 then r1:=FormVariables.my_real_convert(s1,bOk2);  // числовые размеры
   if bOk2 then r2:=FormVariables.my_real_convert(s2,bOk2);  // заданные пользователем
   if bOk2 then r3:=FormVariables.my_real_convert(s3,bOk2);  // с учётом подстановки
   if bOk2 then r4:=FormVariables.my_real_convert(s4,bOk2);  // значений переменных.
   if bOk2 then r5:=FormVariables.my_real_convert(s5,bOk2);
   if bOk2 then r6:=FormVariables.my_real_convert(s6,bOk2);

   if (bOk2) then
   begin
      with Laplas.body[0] do
      begin
         sxS:=s1;  // параметризованные
         syS:=s2;  // геометрические
         sxE:=s3;  // размеры
         syE:=s4;  // заданные
         szS:=s5;  // пользователем
         szE:=s6;

         xS:=r1;  // числовые размеры
         yS:=r2;  // заданные пользователем
         xE:=r3;  // с учётом подстановки
         yE:=r4;  // значений переменных.
         zS:=r5;
         zE:=r6;

         for itekcab23 := 0 to (Laplas.lw-1) do
         begin

         // редактирование стенок геометрические размеры которых зависят
         // от геометрических размеров кабинета.
           case Laplas.wall[itekcab23].cabinet_depend of
              1 : begin
                // minX
                 Laplas.wall[itekcab23].name:='cabinet_minX';
                 Laplas.wall[itekcab23].iPlane:=3;//XY
                 Laplas.wall[itekcab23].xS:=xS;  Laplas.wall[itekcab23].xE:=xS; // числовые значения
                 Laplas.wall[itekcab23].yS:=yS;  Laplas.wall[itekcab23].yE:=yE; // геометрических размеров
                 Laplas.wall[itekcab23].zS:=zS;  Laplas.wall[itekcab23].zE:=zE;
                 Laplas.wall[itekcab23].sxS:=sxS;  Laplas.wall[itekcab23].sxE:=sxS; // параметризованные значения
                 Laplas.wall[itekcab23].syS:=syS;  Laplas.wall[itekcab23].syE:=syE; // геометрических размеров
                 Laplas.wall[itekcab23].szS:=szS;  Laplas.wall[itekcab23].szE:=szE;
                 Laplas.wall[itekcab23].cabinet_depend:=1;
              end;
              2 : begin
                // maxX
                Laplas.wall[itekcab23].name:='cabinet_maxX';
                Laplas.wall[itekcab23].iPlane:=3;//XY
                Laplas.wall[itekcab23].xS:=xE;  Laplas.wall[itekcab23].xE:=xE; // числовые значения
                Laplas.wall[itekcab23].yS:=yS;  Laplas.wall[itekcab23].yE:=yE; // геометрических размеров
                Laplas.wall[itekcab23].zS:=zS;  Laplas.wall[itekcab23].zE:=zE;
                Laplas.wall[itekcab23].sxS:=sxE;  Laplas.wall[itekcab23].sxE:=sxE; // параметризованные значения
                Laplas.wall[itekcab23].syS:=syS;  Laplas.wall[itekcab23].syE:=syE; // геометрических размеров
                Laplas.wall[itekcab23].szS:=szS;  Laplas.wall[itekcab23].szE:=szE;
                Laplas.wall[itekcab23].cabinet_depend:=2;
              end;
              3 : begin
                // minY
                 Laplas.wall[itekcab23].name:='cabinet_minY';
                 Laplas.wall[itekcab23].iPlane:=2;//XZ
                 Laplas.wall[itekcab23].xS:=xS;  Laplas.wall[itekcab23].xE:=xE; // числовые значения
                 Laplas.wall[itekcab23].yS:=yS;  Laplas.wall[itekcab23].yE:=yS; // геометрических размеров
                 Laplas.wall[itekcab23].zS:=zS;  Laplas.wall[itekcab23].zE:=zE;
                 Laplas.wall[itekcab23].sxS:=sxS;  Laplas.wall[itekcab23].sxE:=sxE; // параметризованные значения
                 Laplas.wall[itekcab23].syS:=syS;  Laplas.wall[itekcab23].syE:=syS; // геометрических размеров
                 Laplas.wall[itekcab23].szS:=szS;  Laplas.wall[itekcab23].szE:=szE;
                 Laplas.wall[itekcab23].cabinet_depend:=3;
              end;
              4 : begin
                // maxY
                 Laplas.wall[itekcab23].name:='cabinet_maxY';
                 Laplas.wall[itekcab23].iPlane:=2;//XZ
                 Laplas.wall[itekcab23].xS:=xS;  Laplas.wall[itekcab23].xE:=xE; // числовые значения
                 Laplas.wall[itekcab23].yS:=yE;  Laplas.wall[itekcab23].yE:=yE; // геометрических размеров
                 Laplas.wall[itekcab23].zS:=zS;  Laplas.wall[itekcab23].zE:=zE;
                 Laplas.wall[itekcab23].sxS:=sxS;  Laplas.wall[itekcab23].sxE:=sxE; // параметризованные значения
                 Laplas.wall[itekcab23].syS:=syE;  Laplas.wall[itekcab23].syE:=syE; // геометрических размеров
                 Laplas.wall[itekcab23].szS:=szS;  Laplas.wall[itekcab23].szE:=szE;
                 Laplas.wall[itekcab23].cabinet_depend:=4;
              end;
              5 : begin
                // minZ
                Laplas.wall[itekcab23].name:='cabinet_minZ';
                Laplas.wall[itekcab23].iPlane:=1;//XY
                Laplas.wall[itekcab23].xS:=xS;  Laplas.wall[itekcab23].xE:=xE; // числовые значения
                Laplas.wall[itekcab23].yS:=yS;  Laplas.wall[itekcab23].yE:=yE; // геометрических размеров
                Laplas.wall[itekcab23].zS:=zS;  Laplas.wall[itekcab23].zE:=zS;
                Laplas.wall[itekcab23].sxS:=sxS;  Laplas.wall[itekcab23].sxE:=sxE; // параметризованные значения
                Laplas.wall[itekcab23].syS:=syS;  Laplas.wall[itekcab23].syE:=syE; // геометрических размеров
                Laplas.wall[itekcab23].szS:=szS;  Laplas.wall[itekcab23].szE:=szS;
                Laplas.wall[itekcab23].cabinet_depend:=5;
              end;
              6 : begin
                // maxZ
                Laplas.wall[itekcab23].name:='cabinet_maxZ';
                Laplas.wall[itekcab23].iPlane:=1;//XY
                Laplas.wall[itekcab23].xS:=xS;  Laplas.wall[itekcab23].xE:=xE; // числовые значения
                Laplas.wall[itekcab23].yS:=yS;  Laplas.wall[itekcab23].yE:=yE; // геометрических размеров
                Laplas.wall[itekcab23].zS:=zE;  Laplas.wall[itekcab23].zE:=zE;
                Laplas.wall[itekcab23].sxS:=sxS;  Laplas.wall[itekcab23].sxE:=sxE; // параметризованные значения
                Laplas.wall[itekcab23].syS:=syS;  Laplas.wall[itekcab23].syE:=syE; // геометрических размеров
                Laplas.wall[itekcab23].szS:=szE;  Laplas.wall[itekcab23].szE:=szE;
                Laplas.wall[itekcab23].cabinet_depend:=6;
              end;
           end; // case
         end;
      end;
   end;

   Laplas.ReadyPaint;

   if (not(bOk2)) then
   begin
      bOk:=false;
   end;


    if (bOK) then
    begin
       Close;
    end;

end;

end.
